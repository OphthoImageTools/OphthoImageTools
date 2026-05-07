// ============================================================
//  OCTA Vascular Density & FAZ Analyzer  v3.0
//  Compatible con ImageJ/Fiji
//
//  CAMBIOS v3.1
//  - FAZ_MAX_AREA_MM2: 1.5 → 4.0 mm² (cubre FAZ isquemica)
//  - FAZ_RADIO_BUSQUEDA_UM: 800 → 1500 µm (evita recorte de FAZ grande)
//
//  CAMBIOS v3.0
//  - Preprocesamiento reescrito: sustituido Top-Hat de MorphoLibJ
//    por Subtract Background (nativo ImageJ, sin dependencias)
//  - Corregido bug de densidad 100%: se eliminó imageCalculator 32-bit
//  - Eliminado Fill Holes del paso de vasos (llenaba la FAZ)
//  - FAZ: método primario = Level Sets (Active Contours) como en tu
//    macro original; fallback automático = umbral adaptativo sobre
//    imagen en escala de grises
//
//  TECLAS:
//    F1 → Ejecutar análisis
//    F2 → Configurar parámetros
//    F3 → Exportar tabla a CSV
// ============================================================

// ─────────────────────────────────────────────────────────────
//  PARÁMETROS
// ─────────────────────────────────────────────────────────────
var ESCALA_UM_POR_PIXEL        = 12.0;   // 6 mm / 500 px
var RADIO_CENTRAL_UM           = 500;    // µm – zona perifoveal
var RADIO_PARAFOVEAL_UM        = 1500;   // µm – anillo parafoveal
var FAZ_RADIO_BUSQUEDA_UM      = 1500;   // µm – radio de búsqueda FAZ (aumentado para FAZ isquemica)
var FAZ_MIN_AREA_MM2           = 0.005;
var FAZ_MAX_AREA_MM2           = 4.0;   // aumentado para FAZ isquemica (era 1.5)
var METODO_UMBRAL              = "Otsu";
var ROLLING_BALL_PX            = 20;     // radio background subtraction (px)
var USAR_LEVEL_SETS            = true;   // true = Level Sets (plugin); false = umbral gris
var LS_GREY_THRESHOLD          = 50;     // parámetro grey_value_threshold de Level Sets
var LS_ADVECTION               = 2.20;
var LS_PROPAGATION             = 1;
var LS_CURVATURE               = 1;
var LS_GRAYSCALE               = 30;
var FAZ_UMBRAL_FRACCION        = 0.70;   // para el fallback: % de la media local
var SELECCION_MANUAL_FOVEA     = true;
var GUARDAR_RESULTADOS         = true;
var MODO_LOTE                  = false;

// ─────────────────────────────────────────────────────────────
//  MACRO PRINCIPAL  (F1)
// ─────────────────────────────────────────────────────────────
macro "OCTA Vascular Density & FAZ [F1]" {
    requires("1.53c");
    run("Set Measurements...", "area perimeter shape redirect=None decimal=4");
    setOption("BlackBackground", true);

    if (MODO_LOTE) {
        procesarCarpeta();
    } else {
        if (nImages == 0) {
            ruta = File.openDialog("Seleccionar imagen OCTA");
            open(ruta);
        }
        procesarImagen(getTitle());
    }
    showMessage("Análisis completado", "Resultados en 'OCTA Results'.\nUsa F3 para exportar CSV.");
}

// ─────────────────────────────────────────────────────────────
//  MODO LOTE
// ─────────────────────────────────────────────────────────────
function procesarCarpeta() {
    carpeta  = getDirectory("Seleccionar carpeta con imágenes OCTA");
    if (carpeta == "") return;
    archivos = getFileList(carpeta);
    setBatchMode(true);
    for (i = 0; i < archivos.length; i++) {
        nombre = archivos[i];
        ext = toLowerCase(substring(nombre, lastIndexOf(nombre, ".")));
        if (ext == ".tif" || ext == ".png" || ext == ".jpg" || ext == ".bmp") {
            open(carpeta + nombre);
            procesarImagen(getTitle());
            close();
        }
    }
    setBatchMode(false);
}

// ─────────────────────────────────────────────────────────────
//  PROCESAMIENTO DE UNA IMAGEN
// ─────────────────────────────────────────────────────────────
function procesarImagen(nombreImagen) {

    selectWindow(nombreImagen);
    idOriginal = getImageID();

    // Asegurar 8-bit
    if (bitDepth() != 8) run("8-bit");
    anchoPx = getWidth();
    altoPx  = getHeight();

    // ── PASO 1: Selección manual de la fóvea ─────────────────
    centrX = anchoPx / 2;
    centrY = altoPx  / 2;

    if (SELECCION_MANUAL_FOVEA) {
        selectWindow(nombreImagen);
        setTool("point");
        waitForUser("Selección de fóvea — " + nombreImagen,
            "Haz clic en el CENTRO de la fóvea\n" +
            "(zona avascular central, área más oscura)\n\n" +
            "Pulsa OK para continuar.");
        if (selectionType() == 10) {
            getSelectionCoordinates(xs, ys);
            centrX = xs[0];
            centrY = ys[0];
            run("Select None");
            print("Fóvea: (" + centrX + ", " + centrY + ") px");
        } else {
            print("Aviso: sin selección de punto. Usando centro geométrico.");
        }
    }

    // Radios en píxeles
    rCentral = RADIO_CENTRAL_UM      / ESCALA_UM_POR_PIXEL;
    rPara    = RADIO_PARAFOVEAL_UM   / ESCALA_UM_POR_PIXEL;
    rFAZ     = FAZ_RADIO_BUSQUEDA_UM / ESCALA_UM_POR_PIXEL;

    // ── PASO 2: Preprocesamiento (sin MorphoLibJ) ─────────────
    //
    //  Subtract Background = Rolling Ball = equivalente a Top-Hat
    //  Trabaja en 8-bit in-place, no genera imagen 32-bit intermedia.

    run("Duplicate...", "title=proc_tmp");
    idProc = getImageID();

    // Rolling Ball Background Subtraction (nativo ImageJ)
    run("Subtract Background...", "rolling=" + ROLLING_BALL_PX);

    // CLAHE – realce local de contraste
    run("Enhance Local Contrast (CLAHE)",
        "blocksize=63 histogram=256 maximum=3 mask=*None*");

    // Filtro mediana – eliminar ruido puntual
    run("Median...", "radius=1");

    // ── PASO 3: Segmentación de vasos ────────────────────────
    //
    //  IMPORTANTE: NO usar Fill Holes aquí.
    //  Fill Holes rellena la FAZ (el agujero avascular más grande)
    //  y hace que la densidad suba artificialmente al 100%.

    run("Duplicate...", "title=vasos_mask");
    idVasos = getImageID();

    setAutoThreshold(METODO_UMBRAL + " dark");
    run("Convert to Mask");

    // Solo Opening (erosión + dilatación): elimina artefactos puntales
    // pero NO cierra la FAZ
    run("Open", "stack");

    // ── PASO 4: Densidad vascular ─────────────────────────────

    // 4a) Imagen completa
    selectImage(idVasos);
    getRawStatistics(nTot, mediaTot);
    densTotal = (mediaTot / 255.0) * 100.0;

    // 4b) Zona perifoveal
    makeOval(centrX - rCentral, centrY - rCentral, 2*rCentral, 2*rCentral);
    getRawStatistics(nCentral, mediaCentral);
    densCentral = (mediaCentral / 255.0) * 100.0;
    run("Select None");

    // 4c) Anillo parafoveal
    // Densidad anillo = (vasos en anillo externo - vasos en círculo central) / área anillo
    makeOval(centrX - rPara, centrY - rPara, 2*rPara, 2*rPara);
    getRawStatistics(nParaExt, mediaParaExt);
    run("Select None");

    pixVasosPara  = (mediaParaExt / 255.0) * nParaExt - (mediaCentral / 255.0) * nCentral;
    areaAnillo_px = nParaExt - nCentral;
    densPara      = (pixVasosPara / areaAnillo_px) * 100.0;

    // ── PASO 5: Detección de FAZ ──────────────────────────────
    px2_por_mm2 = 1e6 / (ESCALA_UM_POR_PIXEL * ESCALA_UM_POR_PIXEL);
    fazMinPx2   = FAZ_MIN_AREA_MM2 * px2_por_mm2;
    fazMaxPx2   = FAZ_MAX_AREA_MM2 * px2_por_mm2;

    fazArea_mm2   = 0;
    fazPerim_mm   = 0;
    fazCirc       = 0;
    fazEncontrada = false;
    fazROI        = -1;

    roiManager("reset");
    run("Clear Results");

    if (USAR_LEVEL_SETS) {
        // ── Método A: Level Sets (Active Contours) ────────────
        // Mismo enfoque que tu macro original.
        // Se necesita el plugin "Level Sets" instalado en Fiji.
        // La selección oval inicial se centra en la fóvea con
        // un radio inicial de ~200 µm (ajustable).

        radioInicialPx = 200 / ESCALA_UM_POR_PIXEL; // oval inicial dentro de la FAZ
        selectWindow(nombreImagen);
        makeOval(centrX - radioInicialPx, centrY - radioInicialPx,
                 2*radioInicialPx, 2*radioInicialPx);

        // Ejecutar Level Sets directamente.
        // Si el plugin no está instalado, desactiva USAR_LEVEL_SETS desde F2.
        {
            run("Level Sets",
                "method=[Active Contours] use_level_sets " +
                "grey_value_threshold=" + LS_GREY_THRESHOLD + " " +
                "distance_threshold=0.50 " +
                "advection=" + LS_ADVECTION + " " +
                "propagation=" + LS_PROPAGATION + " " +
                "curvature=" + LS_CURVATURE + " " +
                "grayscale=" + LS_GRAYSCALE + " " +
                "convergence=0.0025 region=outside");

            // "Invert LUT" solo cambia la visualización, no los píxeles.
            // Necesitamos invertir los valores reales para que la FAZ sea
            // blanca (255) y Analyze Particles la encuentre con BlackBackground=true.
            run("Invert");

            // Restringir al radio de búsqueda centrado en la fóvea.
            // Esto excluye zonas de hipoperfusión periféricas fuera del radio.
            makeOval(centrX - rFAZ, centrY - rFAZ, 2*rFAZ, 2*rFAZ);
            run("Make Inverse");
            setColor(0);
            fill();
            run("Select None");

            // Guardar como faz_work para el overlay
            run("Duplicate...", "title=faz_work");
            idFAZwork = getImageID();

            // Analyze Particles: buscar regiones blancas (FAZ) dentro del radio
            // Se elige la región MÁS CERCANA AL CENTRO (punto de fóvea marcado),
            // no la de mayor área, para evitar confusión con hipoperfusiones periféricas.
            run("Analyze Particles...",
                "size=" + fazMinPx2 + "-" + fazMaxPx2 +
                " circularity=0.1-1.0 show=Nothing add");

            if (roiManager("count") > 0) {
                // FAZ = región más cercana al punto de fóvea seleccionado
                minDist = 1e9;
                for (r = 0; r < roiManager("count"); r++) {
                    roiManager("Select", r);
                    getBoundingRect(rx, ry, rw, rh);
                    rcx = rx + rw / 2;
                    rcy = ry + rh / 2;
                    dist = sqrt((rcx - centrX)*(rcx - centrX) + (rcy - centrY)*(rcy - centrY));
                    if (dist < minDist) { minDist = dist; fazROI = r; }
                }
                roiManager("Select", fazROI);
                run("Measure");
                nR = nResults - 1;
                fazArea_mm2   = getResult("Area",   nR) / px2_por_mm2;
                fazPerim_mm   = (getResult("Perim.", nR) * ESCALA_UM_POR_PIXEL) / 1000.0;
                fazCirc       = getResult("Circ.",  nR);
                fazEncontrada = true;
                print("FAZ (Level Sets): " + d2s(fazArea_mm2,3) + " mm²  dist. al centro: " + d2s(minDist,1) + " px");
            } else {
                print("Level Sets no encontró FAZ. Usando método de respaldo...");
                fazMetodoFallback(centrX, centrY, rFAZ,
                                  fazMinPx2, fazMaxPx2, px2_por_mm2);
            }

            // Cerrar la imagen de Level Sets original si sigue abierta
            if (isOpen("Level Sets")) { selectWindow("Level Sets"); close(); }

        }

    } else {
        // ── Método B: Umbral adaptativo sobre imagen en gris ──
        fazMetodoFallback(centrX, centrY, rFAZ,
                          fazMinPx2, fazMaxPx2, px2_por_mm2);
    }

    // ── PASO 6: Guardar resultados ─────────────────────────────
    if (!isOpen("OCTA Results")) Table.create("OCTA Results");
    fila = Table.size("OCTA Results");
    Table.set("Imagen",                  fila, nombreImagen,        "OCTA Results");
    Table.set("Fóvea X (px)",            fila, centrX,              "OCTA Results");
    Table.set("Fóvea Y (px)",            fila, centrY,              "OCTA Results");
    Table.set("Escala (µm/px)",          fila, ESCALA_UM_POR_PIXEL, "OCTA Results");
    Table.set("Densidad Total (%)",      fila, densTotal,           "OCTA Results");
    Table.set("Densidad Central (%)",    fila, densCentral,         "OCTA Results");
    Table.set("Densidad Parafoveal (%)", fila, densPara,            "OCTA Results");
    if (USAR_LEVEL_SETS) fazMetodo = "Level Sets";
    else fazMetodo = "Umbral gris";
    Table.set("FAZ Método",              fila, fazMetodo,            "OCTA Results");
    Table.set("FAZ Encontrada",          fila, fazEncontrada,       "OCTA Results");
    Table.set("FAZ Área (mm²)",          fila, fazArea_mm2,         "OCTA Results");
    Table.set("FAZ Perímetro (mm)",      fila, fazPerim_mm,         "OCTA Results");
    Table.set("FAZ Circularidad",        fila, fazCirc,             "OCTA Results");
    Table.update("OCTA Results");

    // ── PASO 7: Overlay y guardado ────────────────────────────
    if (GUARDAR_RESULTADOS) {
        rutaBase   = getDirectory("image");
        if (rutaBase == "") rutaBase = getDirectory("home");
        nombreBase = replace(replace(replace(
                     nombreImagen, ".tif",""), ".png",""), ".jpg","");

        // Máscara de vasos
        selectImage(idVasos);
        saveAs("PNG", rutaBase + nombreBase + "_vasos.png");

        // Máscara FAZ
        if (isOpen("faz_work")) {
            selectWindow("faz_work");
            saveAs("PNG", rutaBase + nombreBase + "_FAZ.png");
        }

        // Imagen overlay
        selectImage(idOriginal);
        run("Duplicate...", "title=overlay_tmp");
        run("RGB Color");

        setLineWidth(2);
        setColor(255, 220, 0);   // amarillo – zona central
        drawOval(centrX - rCentral, centrY - rCentral, 2*rCentral, 2*rCentral);
        setColor(0, 220, 255);   // cian – parafoveal
        drawOval(centrX - rPara,    centrY - rPara,    2*rPara,    2*rPara);

        if (fazEncontrada && fazROI >= 0) {
            roiManager("Select", fazROI);
            setColor(255, 80, 80);   // rojo – contorno FAZ
            run("Draw", "slice");
        }

        setColor(0, 255, 0);     // verde – punto fóvea
        fillOval(centrX - 5, centrY - 5, 10, 10);

        saveAs("PNG", rutaBase + nombreBase + "_overlay.png");
        selectWindow(nombreBase + "_overlay.png"); close();
    }

    // Limpiar
    roiManager("reset");
    run("Clear Results");
    if (isOpen("proc_tmp"))   { selectWindow("proc_tmp");   close(); }
    if (isOpen("faz_work"))   { selectWindow("faz_work");   close(); }
    if (isOpen("faz_detect")) { selectWindow("faz_detect"); close(); }
    selectImage(idVasos); close();
}

// ─────────────────────────────────────────────────────────────
//  MÉTODO FAZ DE RESPALDO: umbral adaptativo sobre gris
// ─────────────────────────────────────────────────────────────
function fazMetodoFallback(cx, cy, rBusqueda, minPx2, maxPx2, px2mm2) {

    // Trabajar sobre la imagen preprocesada en gris (proc_tmp)
    if (!isOpen("proc_tmp")) {
        // Si ya fue cerrada, usar la original
        selectImage(idOriginal);
        run("Duplicate...", "title=proc_tmp");
        run("Subtract Background...", "rolling=" + ROLLING_BALL_PX);
    } else {
        selectWindow("proc_tmp");
    }

    run("Duplicate...", "title=faz_detect");
    idFAZdetect = getImageID();

    // Fuera del radio de búsqueda → poner blanco (no confundir con FAZ)
    makeOval(cx - rBusqueda, cy - rBusqueda, 2*rBusqueda, 2*rBusqueda);
    run("Make Inverse");
    setColor(255);
    fill();
    run("Select None");

    // Intensidad media dentro del área de búsqueda
    makeOval(cx - rBusqueda, cy - rBusqueda, 2*rBusqueda, 2*rBusqueda);
    getStatistics(aS, meanS);
    run("Select None");

    // La FAZ es más oscura que la media del tejido circundante
    umbralFAZ = meanS * FAZ_UMBRAL_FRACCION;
    setThreshold(0, umbralFAZ);
    run("Convert to Mask");          // FAZ oscura → blanco, resto → negro

    // Limpiar artefactos pequeños
    run("Open", "stack");

    rename("faz_work");
    idFAZwork = getImageID();

    run("Analyze Particles...",
        "size=" + minPx2 + "-" + maxPx2 +
        " circularity=0.1-1.0 show=Nothing add");

    nROIs = roiManager("count");
    if (nROIs == 0) {
        print("AVISO: Método de respaldo tampoco encontró FAZ.");
        print("  Sugerencias:");
        print("  - Verifica que la fóvea esté bien marcada.");
        pr	nt("  - Aumenta FAZ_RADIO_BUSQUEDA_UM o reduce FAZ_UMBRAL_FRACCION (F2).");
        return;
    }

    // La FAZ = región de mayor área dentro del radio de búsqueda
    maxA = 0;
    for (r = 0; r < nROIs; r++) {
        roiManager("Select", r);
        getStatistics(aR);
        if (aR > maxA) { maxA = aR; fazROI = r; }
    }

    roiManager("Select", fazROI);
    run("Measure");
    nR = nResults - 1;
    fazArea_mm2   = getResult("Area",   nR) / px2mm2;
    fazPerim_mm   = (getResult("Perim.", nR) * ESCALA_UM_POR_PIXEL) / 1000.0;
    fazCirc       = getResult("Circ.",  nR);
    fazEncontrada = true;
    print("FAZ (umbral gris, fracción=" + FAZ_UMBRAL_FRACCION + "): " +
          d2s(fazArea_mm2,3) + " mm²");
}



// ─────────────────────────────────────────────────────────────
//  MACRO: Configurar parámetros  (F2)
// ─────────────────────────────────────────────────────────────
macro "Configurar parámetros OCTA [F2]" {
    Dialog.create("Configuración OCTA Analyzer v3");

    Dialog.addMessage("── Escala ─────────────────────────────────────");
    Dialog.addNumber("Escala (µm/pixel):", ESCALA_UM_POR_PIXEL, 2, 7, "µm/px");
    Dialog.addMessage("   6 mm / 500 px = 12.00 µm/px");

    Dialog.addMessage(" ");
    Dialog.addMessage("── Zonas de análisis ──────────────────────────");
    Dialog.addNumber("Radio zona central (perifoveal):", RADIO_CENTRAL_UM,   0, 6, "µm");
    Dialog.addNumber("Radio zona parafoveal:",           RADIO_PARAFOVEAL_UM, 0, 6, "µm");

    Dialog.addMessage(" ");
    Dialog.addMessage("── Preprocesamiento ───────────────────────────");
    Dialog.addNumber("Rolling Ball (background subtraction):", ROLLING_BALL_PX, 0, 4, "px");

    Dialog.addMessage(" ");
    Dialog.addMessage("── Vasos ──────────────────────────────────────");
    items = newArray("Otsu", "Triangle", "MaxEntropy", "Mean", "Yen", "Li");
    Dialog.addChoice("Método umbralización:", items, METODO_UMBRAL);

    Dialog.addMessage(" ");
    Dialog.addMessage("── FAZ ────────────────────────────────────────");
    Dialog.addNumber("Radio búsqueda FAZ (normal ~800, isquemia ~1500):", FAZ_RADIO_BUSQUEDA_UM, 0, 6, "µm");
    Dialog.addNumber("Área mínima FAZ:",    FAZ_MIN_AREA_MM2,      3, 6, "mm²");
    Dialog.addNumber("Área máxima FAZ (normal ~1.5, isquemia ~4.0):", FAZ_MAX_AREA_MM2, 2, 6, "mm²");
    Dialog.addCheckbox("Usar Level Sets (requiere plugin en Fiji)", USAR_LEVEL_SETS);
    Dialog.addMessage("   Parámetros Level Sets:");
    Dialog.addNumber("   grey_value_threshold:", LS_GREY_THRESHOLD, 0, 4, "");
    Dialog.addNumber("   advection:",            LS_ADVECTION,      2, 5, "");
    Dialog.addNumber("   grayscale:",            LS_GRAYSCALE,      0, 4, "");
    Dialog.addMessage("   Parámetro método de respaldo:");
    Dialog.addNumber("   Fracción umbral (0-1):", FAZ_UMBRAL_FRACCION, 2, 5, "");

    Dialog.addMessage(" ");
    Dialog.addMessage("── Opciones ───────────────────────────────────");
    Dialog.addCheckbox("Selección manual de la fóvea", SELECCION_MANUAL_FOVEA);
    Dialog.addCheckbox("Guardar imágenes de resultado", GUARDAR_RESULTADOS);
    Dialog.addCheckbox("Modo lote (procesar carpeta)",  MODO_LOTE);

    Dialog.show();

    ESCALA_UM_POR_PIXEL    = Dialog.getNumber();
    RADIO_CENTRAL_UM       = Dialog.getNumber();
    RADIO_PARAFOVEAL_UM    = Dialog.getNumber();
    ROLLING_BALL_PX        = Dialog.getNumber();
    METODO_UMBRAL          = Dialog.getChoice();
    FAZ_RADIO_BUSQUEDA_UM  = Dialog.getNumber();
    FAZ_MIN_AREA_MM2       = Dialog.getNumber();
    FAZ_MAX_AREA_MM2       = Dialog.getNumber();
    USAR_LEVEL_SETS        = Dialog.getCheckbox();
    LS_GREY_THRESHOLD      = Dialog.getNumber();
    LS_ADVECTION           = Dialog.getNumber();
    LS_GRAYSCALE           = Dialog.getNumber();
    FAZ_UMBRAL_FRACCION    = Dialog.getNumber();
    SELECCION_MANUAL_FOVEA = Dialog.getCheckbox();
    GUARDAR_RESULTADOS     = Dialog.getCheckbox();
    MODO_LOTE              = Dialog.getCheckbox();

    if (USAR_LEVEL_SETS) metodoFAZtxt = "Level Sets";
    else metodoFAZtxt = "Umbral gris";
    showMessage("Configuración guardada",
        "Escala: "         + ESCALA_UM_POR_PIXEL   + " µm/px\n" +
        "Rolling Ball: "   + ROLLING_BALL_PX        + " px\n" +
        "Umbralización: "  + METODO_UMBRAL           + "\n" +
        "Método FAZ: "     + metodoFAZtxt            + "\n" +
        "Fóvea manual: "   + SELECCION_MANUAL_FOVEA);
}

// ─────────────────────────────────────────────────────────────
//  MACRO: Exportar tabla a CSV  (F3)
// ─────────────────────────────────────────────────────────────
macro "Exportar resultados CSV [F3]" {
    if (!isOpen("OCTA Results")) {
        showMessage("Sin resultados", "Ejecuta primero el análisis (F1)."); return;
    }
    ruta = File.saveDialog("Guardar resultados como CSV");
    if (ruta != "") {
        if (!endsWith(ruta, ".csv")) ruta = ruta + ".csv";
        Table.save(ruta, "OCTA Results");
        showMessage("Exportado", "Guardado en:\n" + ruta);
    }
}
