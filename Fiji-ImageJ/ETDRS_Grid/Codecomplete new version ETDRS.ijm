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

//Dibuja circulos y crea ROI
var count = 0;  // Variable to know which circle we are working on

// If you install this macro, it will be mapped to the F2 Key
macro "Draw Concentric Quadrants [F2]" { 
	
	// Clean previous ROIs
	roiManager("Reset");

	// Ask for user to draw a line, to define extents
	makeLine(73, 427, 427, 73)

	// Get the line end points
	getSelectionCoordinates(x,y);

	// Pick up user input
	nCircles   = 1;
	nQuadrants = 4;
	
	// Get angle offset
	offset = atan2(y[1] - y[0], x[1] - x[0]);
	
	// Get Center
	cx = (x[0] + x[1]) /2;
	cy = (y[0] + y[1]) /2;
	
	// Get edge distance
	r = sqrt(pow(x[0] - x[1],2) + pow(y[0] - y[1],2)) / 2;
	
	// Set counter to 0
	count = 0;

	// Make as many concenrtric donuts with quadrants as needed
	for (c=0; c<nCircles; c++) {
		// define inner and outer diameters
		inner = c*(r/nCircles);
		outer = (c+1)*(r/nCircles);

		// Make the quadrants
		doDonut(cx, cy, inner, outer, nQuadrants, offset);
		
	}
}

// This function makes the concentric circles with quadrants
function doDonut(cx, cy, inner, outer, quadrants, offset) {
	count++;
	for (a = 0; a<quadrants; a++) {
		endA = 2*PI / quadrants * a + offset;
		startA   = 2*PI / quadrants * (a+1) + offset;
		step = 30;
		// Build the figure
		pList = newArray(0);
		// Inner circle
		tmp = makeCircle(cx, cy, inner,startA,endA,step);
		pList = Array.concat(pList, tmp);


		// Outer Circle, in the other direction
		tmp = makeCircle(cx, cy, outer, endA, startA,step);
		pList = Array.concat(pList, tmp);
	
		makePolygon(pList);
		Roi.setName("R"+count+"Q"+(a+1));
		roiManager("Add");
	}
}

// Makes a circle from start to end
function makeCircle(cx,cy, radius, startA, endA, nSteps) {
	// We save the coordinates as a 1D array as the macro language does not handle 2D arrays..
	p = newArray((nSteps+1)*2);
	increment = (endA - startA) / nSteps;

	for(i = 0; i <= nSteps; i++) {
		p[2*i]   = radius*cos(startA+i*increment)+cx;
		p[2*i+1] = radius*sin(startA+i*increment)+cy;
	}
	return p;
}

// Convenience function to convert the output array from before into a ROI
function makePolygon(pList) {
	l = pList.length;
	x = newArray(l/2);
	y = newArray(l/2);
	
	for (i=0; i<l/2; i++) {
		x[i] = pList[2*i];
		y[i] = pList[2*i+1];
	}
	
	makeSelection("polygon", x,y);
	
}



