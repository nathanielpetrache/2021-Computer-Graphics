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

int main(int argc, char *argv[]) {
	//
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
	}
	//
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		drawGreyScale(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
