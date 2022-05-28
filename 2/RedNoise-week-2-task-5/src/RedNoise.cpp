#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	float stepSize = (to - from) / (numberOfValues - 1);
	std::vector<float> result;
	for (int step = 0; step < numberOfValues; step++) {
		float current = from + step * stepSize;
		result.push_back(current);
	}	
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues){
	std::vector<float> one = interpolateSingleFloats(from[0], to[0], numberOfValues);
	std::vector<float> two = interpolateSingleFloats(from[1], to[1], numberOfValues);
	std::vector<float> thr = interpolateSingleFloats(from[2], to[2], numberOfValues);
	std::vector<glm::vec3> resultm;
	for (int count = 0; count < numberOfValues; count++) {
		glm::vec3 current = glm::vec3(one[count], two[count], thr[count]);
		resultm.push_back(current);
	}
	return resultm;
	
} 

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

void drawGreyScale(DrawingWindow &window){
	window.clearPixels();
	std::vector<float> values = interpolateSingleFloats(255, 0, window.width);	
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float rgbValues = values[x];
			//std::cout << values[x] << std::endl;
			uint32_t colour = (255 << 24) + (int(rgbValues) << 16) + (int(rgbValues) << 8) + int(rgbValues);
			window.setPixelColour(x, y, colour);
		}
	}
	

}

void drawRainbow(DrawingWindow &window) {

	//corner pixel colours as vec3s
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow
	
	// convert vectors into integers
	
	uint32_t topLeftColour = (255 << 24) + (int(topLeft[0]) << 16) + (int(topLeft[1]) << 8) + int(topLeft[2]);
	uint32_t topRightColour = (255 << 24) + (int(topRight[0]) << 16) + (int(topRight[1]) << 8) + int(topRight[2]);
	uint32_t bottomRightColour = (255 << 24) + (int(bottomRight[0]) << 16) + (int(bottomRight[1]) << 8) + int(bottomRight[2]);
	uint32_t bottomLeftColour = (255 << 24) + (int(bottomLeft[0]) << 16) + (int(bottomLeft[1]) << 8) + int(bottomLeft[2]);
	
	// draw corner pixels	
	// top left
	window.setPixelColour(0, 0, topLeftColour);
	window.setPixelColour(window.width, 0, topRightColour);
	window.setPixelColour(window.width, window.height, bottomRightColour);
	window.setPixelColour(0, window.height, bottomLeftColour);
	
	// interpolate between top and bottom pixels for leftmost column and rightmost column
	//leftmost column
	int x = 0;
	std::vector<glm::vec3> leftColumnColours = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
	for (size_t y = 0; y < window.height; y++) {
		uint32_t currentColour = (255 << 24) + (int(leftColumnColours[y][0]) << 16) + (int(leftColumnColours[y][1]) << 8) + int(leftColumnColours[y][2]);
		window.setPixelColour(x,y, currentColour);
	}
	//rightmost column
	x = window.width-1;
	std::vector<glm::vec3> rightColumnColours = interpolateThreeElementValues(topRight, bottomRight, window.height);
	for (size_t y = 0; y < window.height; y++) {
		uint32_t currentColour = (255 << 24) + (int(rightColumnColours[y][0]) << 16) + (int(rightColumnColours[y][1]) << 8) + int(rightColumnColours[y][2]);
		window.setPixelColour(x,y, currentColour);
	}
	
	// interpolate from left to right for each row
	for (size_t y = 0; y < window.height; y++ ) {
		//get leftmost and rightmost pixels in row
		uint32_t leftmostColour = window.getPixelColour(0, y);
		uint32_t rightmostColour = window.getPixelColour(window.width-1, y);
		/*
		std::cout << "leftmost colour " << leftmostColour;
		std::cout << "rightmost colour " << rightmostColour;
		*/
		//convert colours to glm::vec3 type
		glm::vec3 leftmostColour_v;
		glm::vec3 rightmostColour_v;
		uint32_t redmask = 0x00FF0000; 
		uint32_t gremask = 0x0000FF00;
		uint32_t blumask = 0x000000FF;
		leftmostColour_v[0] = (leftmostColour & redmask) >> 16;
		leftmostColour_v[1] = (leftmostColour & gremask) >> 8;
		leftmostColour_v[2] = (leftmostColour & blumask);
		/*
		std::cout << "leftmost colour red" << leftmostColour_v[0];
		std::cout << "leftmost colour green" << leftmostColour_v[1];
		std::cout << "leftmost colour blue " << leftmostColour_v[2];
		*/
	
		rightmostColour_v[0] = (rightmostColour & redmask) >> 16;
		rightmostColour_v[1] = (rightmostColour & gremask) >> 8;
		rightmostColour_v[2] = (rightmostColour & blumask);
		//interpolate between start and end colours to get vector of vectors of colours 
		std::vector<glm::vec3> currentRowColours = interpolateThreeElementValues(leftmostColour_v, rightmostColour_v, window.width);
		//set each pixel in row to corresponding value in colour vector
		for (size_t x = 0; x < window.width; x++) {
			uint32_t currentColour = (255 << 24) + (int(currentRowColours[x][0]) << 16) + (int(currentRowColours[x][1]) << 8) + int(currentRowColours[x][2]);
		window.setPixelColour(x,y, currentColour);
		}
	}
	
	
}

int main(int argc, char *argv[]) {
	/*
	std::vector<float> result;
	result = interpolateSingleFloats(2.2, 8.5, 7);
	for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	std::cout << std::endl;
	//
	glm::vec3 from(1, 4, 9.2);
	glm::vec3 to(4, 1, 9.8);
	std::vector<glm::vec3> resultm = interpolateThreeElementValues(from, to, 4);
	for (size_t v_i = 0; v_i < resultm.size(); v_i++) {
		for (size_t i=0; i < 3; i++) {
			std::cout << resultm[v_i][i] << " ";
		}
		std::cout << std::endl;
	}*/
	
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		drawRainbow(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
