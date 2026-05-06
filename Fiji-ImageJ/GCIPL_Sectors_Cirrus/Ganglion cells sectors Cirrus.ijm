// ===============================
// 6 Elliptical sectors (Cirrus-like GCIPL)
// ===============================

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
run("Threshold...");
setThreshold(21, 255);

// Run "Analyze Particles" with the specified options
run("Analyze Particles...", "size=4-35 circularity=0.5-1.0 show=[Masks] summarize add");

//Rename
rename("MLCs");
selectImage("MLCs");

roiManager("Reset");
run("Clear Results");

// --- Select foveal center ---
waitForUser("Select foveal center with Point Tool and press OK");
getSelectionCoordinates(xp, yp);
if (xp.length != 1)
   exit("Select exactly one point");

cx = xp[0];
cy = yp[0];

// ellipse
A = 200;   // semi-major axis (x)
B = 167;   // semi-minor axis (y)

// Inner (foveal) ellipse: 100 x 83 px
Ai = 50;   // semi-major
Bi = 41.5; // semi-minor

// Ask for user to draw a line, to define extents
makeLine(cx-A, cy, cx+A, cy);

// Get the line end points
getSelectionCoordinates(x,y);
offset = atan2(y[1] - y[0], x[1] - x[0]);

nSectors = 6;
step = 40;
	
for (s=0; s<nSectors; s++) {
    startA = s*(2*PI/nSectors) + offset;
    endA   = (s+1)*(2*PI/nSectors) + offset;

    p = newArray(0);

    // outer ellipse arc
    for (i=0; i<=step; i++) {
        ang = startA + i*(endA-startA)/step;
        px = cx + A*cos(ang);
        py = cy + B*sin(ang);
        p = Array.concat(p, newArray(px,py));
    }

  // --- inner ellipse arc (reverse direction) ---
for (i=step; i>=0; i--) {
    ang = startA + i*(endA-startA)/step;
    px = cx + Ai*cos(ang);
    py = cy + Bi*sin(ang);
    p = Array.concat(p, newArray(px,py));
}

    makePolygonFromArray(p);
    Roi.setName("GCIPL_S"+(s+1));
    roiManager("Add");
}

// ---------- COUNT PARTICLES PER SECTOR ----------
for (i=0; i<roiManager("count"); i++) {
    roiManager("Select", i);
    run("Analyze Particles...", "size=4-35 circularity=0.5-1.0 clear summarize");
}

// ---------- FUNCTIONS ----------
function makePolygonFromArray(p) {
    n = p.length/2;
    x = newArray(n);
    y = newArray(n);
    for (i=0; i<n; i++) {
        x[i] = p[2*i];
        y[i] = p[2*i+1];
    }
    makeSelection("polygon", x,y);
}


