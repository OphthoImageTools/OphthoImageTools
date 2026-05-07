// Macro para dibujar una elipse con uno de los focos en un punto seleccionado
run("Clear Results");
waitForUser("Select the center of the optic nerve (clic on point tool) and press OK.");

// Obtener coordenadas del punto seleccionado
getSelectionCoordinates(xpoints, ypoints);
if (xpoints.length != 1) {
    exit("Please, select only one point.");
}
fx = xpoints[0];
fy = ypoints[0];

// Solicitar tamaño elipse
a = getNumber("Introduce width", 100);
b = getNumber("introduce height", 50);

// Dibujar la elipse
makeOval(fx-(a/2),fy-(b/2),a,b);

//Recortar la elipse
run("Crop");
run("Duplicate...");
run("8-bit");
run("6 shades");
run("Duplicate...");
run("Threshold...");
setThreshold(80, 255);
setOption("BlackBackground", true);
run("Convert to Mask");
run("Analyze Particles...", "summarize add");
