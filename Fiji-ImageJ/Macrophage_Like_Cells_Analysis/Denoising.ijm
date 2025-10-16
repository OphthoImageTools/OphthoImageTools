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
