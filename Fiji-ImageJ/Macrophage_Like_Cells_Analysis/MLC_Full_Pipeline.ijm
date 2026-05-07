// ============================================================
//  MLC Full Pipeline — COMPLETO
//  v3: preprocesado (registro, DoG, FFT) + analisis (ETDRS, NND)
//      fusionados en un unico macro
//
//  Flujo:
//    1. Dialogo de parametros
//    2. Registro y average OCTA
//    3. Transform y average OCTR
//    4. Filtro DoG (Difference of Gaussians)
//    5. FFT + eliminacion de artefactos de linea
//    6. Inverse FFT + umbralización Triangle estabilizada
//    7. Analyze Particles
//    8. Grid ETDRS (9 regiones, OD/OS)
//    9. Calculo de NND global
//   10. Estadisticas zonales por region ETDRS
//   11. Grafico + guardado de resultados
// ============================================================

macro "MLC Full Pipeline" {

    // ----------------------------------------------------------
    // 1. DIALOGO DE PARAMETROS
    // ----------------------------------------------------------
    Dialog.create("MLC Full Pipeline — Parametros");

    Dialog.addMessage("=== Rutas ===");
    Dialog.addString("Carpeta base (OD o OS):",
        "/Users/estercarrenosalas/Documents/OD/", 60);

    Dialog.addMessage("=== Ojo ===");
    Dialog.addChoice("Ojo:", newArray("OD (ojo derecho)", "OS (ojo izquierdo)"));

    Dialog.addMessage("=== Escala ===");
    Dialog.addNumber("Tamano del campo (mm):", 6);
    Dialog.addNumber("Resolucion (pixeles):", 500);

    Dialog.addMessage("=== Filtro DoG ===");
    Dialog.addNumber("Sigma 1 (blur suave):", 1);
    Dialog.addNumber("Sigma 2 (blur fuerte):", 3);

    Dialog.addMessage("=== Umbral Triangle ===");
    Dialog.addNumber("Blur pre-threshold (sigma, 0 = sin blur):", 1);

    Dialog.addMessage("=== Analyze Particles ===");
    Dialog.addNumber("Tamano minimo (px²):", 4);
    Dialog.addNumber("Tamano maximo (px²):", 35);
    Dialog.addNumber("Circularidad minima (0-1):", 0.5);

    Dialog.show();

    // Recoger valores
    basePath  = Dialog.getString();
    if (!endsWith(basePath, "/")) basePath = basePath + "/";

    eyeChoice = Dialog.getChoice();
    field_mm  = Dialog.getNumber();
    res_px    = Dialog.getNumber();
    sigma1    = Dialog.getNumber();
    sigma2    = Dialog.getNumber();
    blurPre   = Dialog.getNumber();
    sizeMin   = Dialog.getNumber();
    sizeMax   = Dialog.getNumber();
    circMin   = Dialog.getNumber();

    // Rutas derivadas
    srcOCTA      = basePath + "OCTA";
    outOCTA      = basePath + "OutputOCTA";
    srcOCTR      = basePath + "OCTR";
    outOCTR      = basePath + "OutputOCTR";
    transformDir = basePath + "Transform";
    avgDir       = basePath + "Average/";

    // Parametros derivados
    isOD      = startsWith(eyeChoice, "OD");
    if (isOD) { eyeStr = "OD"; } else { eyeStr = "OS"; }
    pxPerMm   = res_px / field_mm;
    um_per_px = (field_mm * 1000.0) / res_px;
    area_mm2  = field_mm * field_mm;

    // ----------------------------------------------------------
    // 2. REGISTRO Y AVERAGE DE OCTA
    // ----------------------------------------------------------
    run("Register Virtual Stack Slices",
        "source=[" + srcOCTA + "] " +
        "output=[" + outOCTA + "] " +
        "feature=Rigid registration=[Elastic              -- bUnwarpJ splines                    ] save");

    selectImage("Registered OCTA");
    run("RGB Color");
    run("AvgNoiseRmvr ");
    saveAs("PNG", avgDir + "averageOCTA.png");
    selectImage("Registered OCTA");
    close();

    // ----------------------------------------------------------
    // 3. TRANSFORM Y AVERAGE DE OCTR
    // ----------------------------------------------------------
    run("Transform Virtual Stack Slices",
        "source=[" + srcOCTR + "] " +
        "output=[" + outOCTR + "] " +
        "transforms=[" + transformDir + "] interpolate");

    selectImage("Registered OCTR");
    run("RGB Color");
    run("AvgNoiseRmvr ");
    saveAs("PNG", avgDir + "averageOCTR.png");
    selectImage("Registered OCTR");
    close();

    // ----------------------------------------------------------
    // 4. DOG (Difference of Gaussians)
    // ----------------------------------------------------------
    selectImage("averageOCTR.png");
    run("Duplicate...", "title=blur1");
    run("Gaussian Blur...", "sigma=" + sigma1);

    selectImage("averageOCTR.png");
    run("Gaussian Blur...", "sigma=" + sigma2);

    imageCalculator("Subtract create", "blur1", "averageOCTR.png");
    dogTitle = getTitle();
    selectImage("blur1");
    close();

    // ----------------------------------------------------------
    // 5. FILTRADO FFT
    // ----------------------------------------------------------
    selectImage(dogTitle);
    run("FFT");

    getDimensions(width, height, channels, slices, frames);
    centerX = width  / 2;
    centerY = height / 2;

    setColor(0, 0, 0);
    fillRect(centerX - 3, 0, 7, height);          // linea vertical
    fillRect(0, centerY - 3, width, 7);           // linea horizontal
    fillOval(centerX - 10, centerY - 10, 20, 20); // componente DC

    run("Inverse FFT");
    ifftTitle = getTitle();

    selectImage(dogTitle);
    close();

    saveAs("PNG", avgDir + "averageOCTRclean.png");

    // ----------------------------------------------------------
    // 6. UMBRALIZACIÓN TRIANGLE ESTABILIZADA
    // ----------------------------------------------------------
    selectImage("averageOCTRclean.png");
    run("8-bit");

    // Suavizado previo: reduce sensibilidad de Triangle a 1 punto
    if (blurPre > 0)
        run("Gaussian Blur...", "sigma=" + blurPre);

    setAutoThreshold("Triangle dark no-reset");
    run("Convert to Mask");

    // ----------------------------------------------------------
    // 7. ANALYZE PARTICLES
    // ----------------------------------------------------------
    run("Set Measurements...", "area centroid shape display redirect=None decimal=3");
    run("Clear Results");
    roiManager("Reset");

    run("Analyze Particles...",
        "size=" + sizeMin + "-" + sizeMax +
        " circularity=" + circMin + "-1.00" +
        " show=Masks display summarize exclude clear add");

    rename("MLCs_mask");

    n = nResults;
    if (n == 0)
        exit("No se detectaron particulas.\nRevisa los parametros de umbral o tamano.");

    print("\\Clear");
    print("=== MLC Full Pipeline (" + eyeStr + ") ===");
    print("Imagen original : averageOCTR.png");
    print("Celulas totales : " + n);
    print("Escala          : " + d2s(um_per_px, 2) + " um/px");

    // ----------------------------------------------------------
    // 8. GRID ETDRS — seleccion de fovea
    // ----------------------------------------------------------
    r1 = 0.5 * pxPerMm;
    r2 = 1.5 * pxPerMm;
    r3 = 3.0 * pxPerMm;

    supS = -3*PI/4;  supE = -PI/4;
    infS =    PI/4;  infE =  3*PI/4;
    rgtS =   -PI/4;  rgtE =  PI/4;
    lftS =  3*PI/4;  lftE =  5*PI/4;

    if (isOD) {
        tmpS = rgtS; tmpE = rgtE;
        nasS = lftS; nasE = lftE;
    } else {
        tmpS = lftS; tmpE = lftE;
        nasS = rgtS; nasE = rgtE;
    }

    selectImage("MLCs_mask");
    setTool("point");
    waitForUser("Centro de la fovea",
        "Haz clic en el CENTRO DE LA FOVEA.\nPulsa OK cuando estes listo.");

    getSelectionCoordinates(xpts, ypts);
    if (xpts.length != 1)
        exit("Selecciona exactamente UN punto sobre la fovea.");
    cx = xpts[0];
    cy = ypts[0];

    print("");
    print("Centro fovea : (" + d2s(cx,1) + ", " + d2s(cy,1) + ") px");
    print("Radio r1 (central)  : " + d2s(r1,1) + " px = 0.5 mm");
    print("Radio r2 (interno)  : " + d2s(r2,1) + " px = 1.5 mm");
    print("Radio r3 (externo)  : " + d2s(r3,1) + " px = 3.0 mm");

    roiManager("Reset");

    makeOval(cx - r1, cy - r1, 2*r1, 2*r1);
    Roi.setName("C");
    roiManager("Add");

    addSector("IS", cx, cy, r1, r2, supS, supE);
    addSector("IT", cx, cy, r1, r2, tmpS, tmpE);
    addSector("II", cx, cy, r1, r2, infS, infE);
    addSector("IN", cx, cy, r1, r2, nasS, nasE);

    addSector("OS", cx, cy, r2, r3, supS, supE);
    addSector("OT", cx, cy, r2, r3, tmpS, tmpE);
    addSector("OI", cx, cy, r2, r3, infS, infE);
    addSector("ON", cx, cy, r2, r3, nasS, nasE);

    roiManager("Show All");
    run("Select None");

    // ----------------------------------------------------------
    // 9. CALCULO DE NND GLOBAL
    // ----------------------------------------------------------
    xc = newArray(n);
    yc = newArray(n);
    for (i = 0; i < n; i++) {
        xc[i] = getResult("X", i);
        yc[i] = getResult("Y", i);
    }

    nnd_um = newArray(n);
    nn_idx = newArray(n);

    for (i = 0; i < n; i++) {
        minD = 1e15;
        best = -1;
        for (j = 0; j < n; j++) {
            if (i == j) continue;
            dx = xc[i] - xc[j];
            dy = yc[i] - yc[j];
            d  = sqrt(dx*dx + dy*dy);
            if (d < minD) { minD = d; best = j; }
        }
        nnd_um[i] = minD * um_per_px;
        nn_idx[i] = best;
    }

    for (i = 0; i < n; i++) {
        setResult("NND_um",   i, nnd_um[i]);
        setResult("NN_index", i, nn_idx[i]);
    }
    updateResults();

    // ----------------------------------------------------------
    // 10. ASIGNAR CELULAS A REGIONES ETDRS
    // ----------------------------------------------------------
    rNames  = newArray("C",  "IS",  "IT",  "II",  "IN",  "OS",  "OT",  "OI",  "ON");
    rLabels = newArray("Central",
                       "Inner Sup", "Inner Temp", "Inner Inf", "Inner Nas",
                       "Outer Sup", "Outer Temp", "Outer Inf", "Outer Nas");

    areaC  = PI * 0.5 * 0.5;
    areaI  = (PI * 1.5 * 1.5 - PI * 0.5 * 0.5) / 4.0;
    areaO  = (PI * 3.0 * 3.0 - PI * 1.5 * 1.5) / 4.0;
    rAreas = newArray(areaC, areaI, areaI, areaI, areaI,
                             areaO, areaO, areaO, areaO);
    nR = rNames.length;

    nROIs = roiManager("count");
    roiForRegion = newArray(nR);
    for (r = 0; r < nR; r++) roiForRegion[r] = -1;
    for (ri = 0; ri < nROIs; ri++) {
        roiManager("select", ri);
        rn = Roi.getName();
        for (r = 0; r < nR; r++) {
            if (rn == rNames[r]) roiForRegion[r] = ri;
        }
    }

    assignment = newArray(n);
    for (i = 0; i < n; i++) assignment[i] = -1;

    for (r = 0; r < nR; r++) {
        if (roiForRegion[r] < 0) continue;
        roiManager("select", roiForRegion[r]);
        for (i = 0; i < n; i++) {
            if (assignment[i] == -1) {
                if (Roi.contains(xc[i], yc[i])) {
                    assignment[i] = r;
                }
            }
        }
    }

    for (i = 0; i < n; i++) {
        r = assignment[i];
        if (r >= 0) { setResult("ETDRS_region", i, rNames[r]); }
        else        { setResult("ETDRS_region", i, "Fuera"); }
    }
    updateResults();

    // ----------------------------------------------------------
    // 11. ESTADISTICAS ZONALES
    // ----------------------------------------------------------
    rCount = newArray(nR);
    rSum   = newArray(nR);
    rMin   = newArray(nR);
    rMax   = newArray(nR);
    for (r = 0; r < nR; r++) { rMin[r] = 1e15; rMax[r] = 0; }

    for (i = 0; i < n; i++) {
        r = assignment[i];
        if (r < 0) continue;
        rCount[r]++;
        rSum[r] += nnd_um[i];
        if (nnd_um[i] < rMin[r]) rMin[r] = nnd_um[i];
        if (nnd_um[i] > rMax[r]) rMax[r] = nnd_um[i];
    }

    rMean = newArray(nR);
    rSD   = newArray(nR);
    rMed  = newArray(nR);
    rDens = newArray(nR);

    for (r = 0; r < nR; r++) {
        if (rCount[r] == 0) continue;
        rMean[r] = rSum[r] / rCount[r];
        rDens[r] = rCount[r] / rAreas[r];

        tmp = newArray(rCount[r]);
        idx = 0;
        for (i = 0; i < n; i++) {
            if (assignment[i] == r) { tmp[idx] = nnd_um[i]; idx++; }
        }

        ssq = 0;
        for (k = 0; k < rCount[r]; k++) {
            d = tmp[k] - rMean[r];
            ssq += d * d;
        }
        if (rCount[r] > 1) rSD[r] = sqrt(ssq / (rCount[r] - 1));

        Array.sort(tmp);
        nk = rCount[r];
        if (nk % 2 == 0) { rMed[r] = (tmp[nk/2 - 1] + tmp[nk/2]) / 2.0; }
        else             { rMed[r] = tmp[floor(nk/2)]; }
    }

    // ----------------------------------------------------------
    // 12. LOG DE RESULTADOS
    // ----------------------------------------------------------
    gSum = 0; gMin = 1e15; gMax = 0;
    for (i = 0; i < n; i++) {
        gSum += nnd_um[i];
        if (nnd_um[i] < gMin) gMin = nnd_um[i];
        if (nnd_um[i] > gMax) gMax = nnd_um[i];
    }
    gMean = gSum / n;
    gSSQ  = 0;
    for (i = 0; i < n; i++) { d = nnd_um[i] - gMean; gSSQ += d*d; }
    gSD   = sqrt(gSSQ / (n - 1));

    sortedG = Array.copy(nnd_um);
    Array.sort(sortedG);
    if (n % 2 == 0) { gMed = (sortedG[n/2-1] + sortedG[n/2]) / 2.0; }
    else            { gMed = sortedG[floor(n/2)]; }

    print("");
    print("--- Estadisticas globales ---");
    print("NND media    : " + d2s(gMean, 1) + " um");
    print("NND mediana  : " + d2s(gMed,  1) + " um");
    print("NND SD       : " + d2s(gSD,   1) + " um");
    print("Densidad     : " + d2s(n / area_mm2, 1) + " cells/mm2");
    print("");
    print("--- Estadisticas por region ETDRS ---");
    print("Region           N      Dens(c/mm2)  NND media(um)  NND SD(um)  NND mediana(um)");
    print("---------------------------------------------------------------------------------");

    for (r = 0; r < nR; r++) {
        if (rCount[r] == 0) {
            print(padR(rLabels[r], 17) + "0      —            —              —           —");
        } else {
            print(
                padR(rLabels[r], 17) +
                padL("" + rCount[r],    5) + "  " +
                padL(d2s(rDens[r], 1), 11) + "  " +
                padL(d2s(rMean[r], 1), 13) + "  " +
                padL(d2s(rSD[r],   1), 10) + "  " +
                d2s(rMed[r], 1)
            );
        }
    }

    nFuera = 0;
    for (i = 0; i < n; i++) if (assignment[i] < 0) nFuera++;
    if (nFuera > 0) print("\nCelulas fuera de la grid ETDRS: " + nFuera);

    print("\nTabla guardada en: " + avgDir);

    // ----------------------------------------------------------
    // 13. HEATMAP DE DENSIDAD CELULAR (15x15 px por celda)
    // ----------------------------------------------------------
    // Deseleccionar ROIs ETDRS: si queda uno activo, Duplicate
    // solo copia esa region en lugar de la imagen completa
    roiManager("Deselect");
    run("Select None");

    selectImage("MLCs_mask");
    getDimensions(width, height, channels, slices, frames);
    run("Duplicate...", "title=Graph");
    run("RGB Color");

    roiSize = 15;
    nROIsX  = floor(width  / roiSize);
    nROIsY  = floor(height / roiSize);

    selectImage("Graph");
    for (roiX = 0; roiX < nROIsX; roiX++) {
        for (roiY = 0; roiY < nROIsY; roiY++) {
            startX = roiX * roiSize;
            startY = roiY * roiSize;

            // Contar centroides de celulas que caen en este ROI
            // usando los arrays xc[], yc[] ya calculados
            cellsInROI = 0;
            for (i = 0; i < n; i++) {
                if (xc[i] >= startX && xc[i] < startX + roiSize &&
                    yc[i] >= startY && yc[i] < startY + roiSize) {
                    cellsInROI++;
                }
            }

            makeRectangle(startX, startY, roiSize, roiSize);

            if (cellsInROI >= 4) {
                setColor(255, 255, 255); // Blanco  — muy alta densidad
            } else if (cellsInROI >= 3) {
                setColor(255, 0,   0  ); // Rojo    — alta densidad
            } else if (cellsInROI >= 2) {
                setColor(255, 255, 0  ); // Amarillo — densidad media
            } else if (cellsInROI >= 1) {
                setColor(0,   255, 0  ); // Verde   — densidad baja
            } else {
                setColor(0,   0,   255); // Azul    — sin celulas
            }
            fill();
        }
    }
    run("Select None");
    updateDisplay();

    // ----------------------------------------------------------
    // 14. GRAFICO NND POR REGION
    // ----------------------------------------------------------
    xPos  = newArray(nR);
    yMean = newArray(nR);
    ySD   = newArray(nR);
    maxY  = 0;
    for (r = 0; r < nR; r++) {
        xPos[r]  = r + 1;
        yMean[r] = rMean[r];
        ySD[r]   = rSD[r];
        if (rMean[r] + rSD[r] > maxY) maxY = rMean[r] + rSD[r];
    }

    Plot.create("NND por region ETDRS — " + eyeStr,
                "Region (1=C, 2-5=Inner, 6-9=Outer)", "NND media (um)");
    Plot.setColor("blue");
    Plot.add("circles", xPos, yMean);
    Plot.setColor("#6baed6");
    Plot.add("error bars", xPos, yMean, ySD);
    Plot.setLimits(0, nR + 1, 0, maxY * 1.3 + 20);
    Plot.show();
    plotTitle = getTitle();

    // ----------------------------------------------------------
    // 15. IMAGEN COMBINED (OCTA + MLCs en rojo)
    // ----------------------------------------------------------
    selectImage("MLCs_mask");
    run("Invert");
    run("Red");
    run("Invert");
    imageCalculator("Add create", "averageOCTA.png", "MLCs_mask");
    combinedTitle = getTitle();

    // ----------------------------------------------------------
    // 16. GUARDAR LOG EN FICHERO TXT
    // ----------------------------------------------------------
    logPath = avgDir + "Log_" + eyeStr + ".txt";
    f = File.open(logPath);
    print(f, "=== MLC Full Pipeline (" + eyeStr + ") ===");
    print(f, "Imagen original : averageOCTR.png");
    print(f, "Celulas totales : " + n);
    print(f, "Escala          : " + d2s(um_per_px, 2) + " um/px");
    print(f, "");
    print(f, "Centro fovea : (" + d2s(cx,1) + ", " + d2s(cy,1) + ") px");
    print(f, "Radio r1 (central)  : " + d2s(r1,1) + " px = 0.5 mm");
    print(f, "Radio r2 (interno)  : " + d2s(r2,1) + " px = 1.5 mm");
    print(f, "Radio r3 (externo)  : " + d2s(r3,1) + " px = 3.0 mm");
    print(f, "");
    print(f, "--- Estadisticas globales ---");
    print(f, "NND media    : " + d2s(gMean, 1) + " um");
    print(f, "NND mediana  : " + d2s(gMed,  1) + " um");
    print(f, "NND SD       : " + d2s(gSD,   1) + " um");
    print(f, "Densidad     : " + d2s(n / area_mm2, 1) + " cells/mm2");
    print(f, "");
    print(f, "--- Estadisticas por region ETDRS ---");
    print(f, "Region           N      Dens(c/mm2)  NND media(um)  NND SD(um)  NND mediana(um)");
    print(f, "---------------------------------------------------------------------------------");
    for (r = 0; r < nR; r++) {
        if (rCount[r] == 0) {
            print(f, padR(rLabels[r], 17) + "0      -            -              -           -");
        } else {
            print(f,
                padR(rLabels[r], 17) +
                padL("" + rCount[r],    5) + "  " +
                padL(d2s(rDens[r], 1), 11) + "  " +
                padL(d2s(rMean[r], 1), 13) + "  " +
                padL(d2s(rSD[r],   1), 10) + "  " +
                d2s(rMed[r], 1)
            );
        }
    }
    if (nFuera > 0) print(f, "\nCelulas fuera de la grid ETDRS: " + nFuera);
    print(f, "\nTabla guardada en: " + avgDir);
    File.close(f);

    // ----------------------------------------------------------
    // 17. GUARDAR IMAGENES Y CERRAR TODO
    // ----------------------------------------------------------
    // Combined
    selectImage(combinedTitle);
    saveAs("PNG", avgDir + "combined.png");
    close();

    // MLCs
    selectImage("MLCs_mask");
    saveAs("PNG", avgDir + "MLCs.png");
    close();

    // Heatmap de densidad
    selectImage("Graph");
    saveAs("PNG", avgDir + "Graph.png");
    close();

    // Grafico NND
    selectImage(plotTitle);
    saveAs("PNG", avgDir + "NND_plot_" + eyeStr + ".png");
    close();

    // CSV con resultados por celula
    saveAs("Results", avgDir + "MLC_NND_ETDRS_" + eyeStr + ".csv");

    // Cerrar el resto de imagenes abiertas
    run("Close All");

    showMessage("Pipeline completado (" + eyeStr + ")",
        "Celulas detectadas : " + n + "\n" +
        "NND media global   : " + d2s(gMean, 1) + " um\n" +
        "Densidad global    : " + d2s(n / area_mm2, 1) + " cells/mm2\n\n" +
        "Archivos guardados en:\n" + avgDir + "\n\n" +
        "  - averageOCTA.png\n" +
        "  - averageOCTR.png\n" +
        "  - averageOCTRclean.png\n" +
        "  - MLCs.png\n" +
        "  - Graph.png\n" +
        "  - combined.png\n" +
        "  - NND_plot_" + eyeStr + ".png\n" +
        "  - MLC_NND_ETDRS_" + eyeStr + ".csv\n" +
        "  - Log_" + eyeStr + ".txt");

} // fin macro


// ============================================================
//  FUNCIONES AUXILIARES
// ============================================================

function addSector(name, cx, cy, inner, outer, startA, endA) {
    nSteps = 60;
    pList  = newArray(0);
    tmp    = arcPoints(cx, cy, inner, startA, endA, nSteps);
    pList  = Array.concat(pList, tmp);
    tmp    = arcPoints(cx, cy, outer, endA, startA, nSteps);
    pList  = Array.concat(pList, tmp);

    np = pList.length / 2;
    xp = newArray(np);
    yp = newArray(np);
    for (i = 0; i < np; i++) {
        xp[i] = pList[2*i];
        yp[i] = pList[2*i+1];
    }
    makeSelection("polygon", xp, yp);
    Roi.setName(name);
    roiManager("Add");
}

function arcPoints(cx, cy, radius, startA, endA, nSteps) {
    p   = newArray((nSteps + 1) * 2);
    inc = (endA - startA) / nSteps;
    for (i = 0; i <= nSteps; i++) {
        p[2*i]   = radius * cos(startA + i * inc) + cx;
        p[2*i+1] = radius * sin(startA + i * inc) + cy;
    }
    return p;
}

function padR(s, w) {
    while (lengthOf(s) < w) s = s + " ";
    return s;
}

function padL(s, w) {
    while (lengthOf(s) < w) s = " " + s;
    return s;
}
