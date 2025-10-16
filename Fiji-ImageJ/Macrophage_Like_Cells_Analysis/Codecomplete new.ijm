// Registrar y average de OCTA
run("Register Virtual Stack Slices", "source=[/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/OCTA] output=[/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/OutputOCTA] feature=Rigid registration=[Elastic              -- bUnwarpJ splines                    ] save");
selectImage("Registered OCTA");
run("RGB Color");
run("AvgNoiseRmvr ");
saveAs("PNG", "/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/Average/averageOCTA.png");
selectImage("Registered OCTA");
close();
selectImage("averageOCTA.png");

//Transform y average de OCTR
run("Transform Virtual Stack Slices", "source=[/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/OCTR] output=[/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/OutputOCTR] transforms=[/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/Transform] interpolate");
selectImage("Registered OCTR");
run("RGB Color");
run("AvgNoiseRmvr ");
saveAs("PNG", "/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/Average/averageOCTR.png");
selectImage("Registered OCTR");
close();
selectImage("averageOCTR.png");

// Perform FFT on the current open image
run("FFT");

// Get image dimensions
getDimensions(width, height, channels, slices, frames);

// Calculate center coordinates
centerX = width / 2;
centerY = height / 2;

// Set cross length
crossLength = Math.min(width, height)/2;

// Draw a black line on the FFT image
drawRect(centerX - 3, 0, 7, width); 
setColor(0, 0, 0);
fillRect(centerX - 3, 0, 7, width);

// Draw second black line on the FFT image
drawRect(0, centerY - 3, height, 7); 
setColor(0, 0, 0);
fillRect(0, centerY - 3, height, 7);

//Draw black circle on the FFT image
drawOval(centerX - 10, centerY - 10, 20, 20);
setColor(0, 0, 0);
fillOval(centerX - 10, centerY - 10, 20, 20);

// Calculate the inverse FFT
run("Inverse FFT");

// Convert the image to 8-bit grayscale
run("8-bit");

// Apply the threshold
setThreshold(21, 255);
run("Apply LUT");

// Run "Analyze Particles" with the specified options
run("Analyze Particles...", "size=4-35 circularity=0.5-1.0 show=[Masks] summarize add");

// Get the results table
results = getResult("Results");

// Set the particle size
particleSize = 4-35;
particleCircularity= 0.5-1.0

// Get the particle count
nParticles = nResults;

//Rename
rename("MLCs");

// Duplicate
run("Duplicate...", "title=Graph")

// Name duplicate as Graph

// Get image dimensions
getDimensions(width, height, channels, slices, frames);

// Define ROI size
roiSize = 15;

// Calculate the number of ROIs in each dimension
nROIsX = width / roiSize;
nROIsY = height / roiSize;

// Convert the image to RGB color mode
run("16_colors");

// Iterate over each ROI
for (roiX = 0; roiX < nROIsX; roiX++) {
    for (roiY = 0; roiY < nROIsY; roiY++) {
        // Calculate ROI coordinates
        startX = roiX * roiSize;
        startY = roiY * roiSize;

        // Set the ROI on the duplicate image
        makeRectangle(startX, startY, roiSize, roiSize);
        roiManager("Add");

        // Run "Analyze Particles" with the specified options on the thresholded image
        setAutoThreshold();
        run("Analyze Particles...", "size=4-35 circularity=0.5-1.0 show summarize add");

        // Get the particle count and set the color for ROI
        v = Table.get("Count");
         if (v >= 4) {
            setColor(255, 255, 255); // White
        } else if (v >= 4) {
            setColor(255, 255, 255); // White
        } else if (v >= 3) {
            setColor(255, 0, 0); // Red
        } else if (v >= 2) {
            setColor(255, 255, 0); // Yellow
        } else if (v >= 1) {
            setColor(0, 255, 0); // Green
        } else {
            setColor(0, 0, 255); // Blue
        }
        run("Fill");       
    }
}

// Update the image display
updateDisplay();

// Cambiar a 16 colores
run("16_colors");

//Select and conver the MLCs
selectImage("MLCs");
run("Invert")
run("Red")
run("Invert")

// Add the two images
imageCalculator("Add create", "averageOCTA.png", "MLCs")

// Guardar imágenes y cerrar el resto
saveAs("PNG", "/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/Average/combined.png");
close();
saveAs("PNG", "/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/Average/MLCs.png");
close();
saveAs("PNG", "/Users/estercarrenosalas/Documents/Documentos - MacBook Air de Ester/Macrofagos/DPV/OS no DPV/Average/Graph.png");
close();
run("Close All");
