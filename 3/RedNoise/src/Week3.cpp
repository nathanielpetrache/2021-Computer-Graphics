#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <CanvasPoint.h>
#include <Colour.h>
#include <time.h>
#include <CanvasTriangle.h>
#include <algorithm>
#include <TextureMap.h>

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

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c) {
	//window.clearPixels();
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(std::abs(xDiff), std::abs(yDiff));
	float xStepSize = xDiff / numberOfSteps;
	float yStepSize = yDiff / numberOfSteps;
	// Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
	for (float step = 0.0; step < numberOfSteps; step++) {
		float x = from.x + (step * xStepSize);
		float y = from.y + (step * yStepSize);
		//Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
		uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + (c.blue);
		window.setPixelColour(std::round(x), std::round(y), colour);
	}	
}

void drawStrokedTriangle(DrawingWindow &window, CanvasTriangle t, Colour c) {
	CanvasPoint i = t.vertices[0];
	CanvasPoint j = t.vertices[1];
	CanvasPoint k = t.vertices[2];
	drawLine(window, i, j, c);
	drawLine(window, i, k, c);
	drawLine(window, j, k, c);
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle t, Colour c) {
	CanvasPoint i = t.vertices[0];
	CanvasPoint j = t.vertices[1];
	CanvasPoint k = t.vertices[2];
	std::vector<CanvasPoint> points;
	//sort points
	if ((i.y <= j.y) && (i.y <= k.y)) {
		points.push_back(i);
		if (j.y <= k.y) {
			points.push_back(j);
			points.push_back(k);
		} else {
			points.push_back(k);
			points.push_back(j);
		}
	} else if ((j.y <= i.y) && (j.y <= k.y)) {
		points.push_back(j);
		if (i.y <= k.y) {
			points.push_back(i);
			points.push_back(k);
		} else {
			points.push_back(k);
			points.push_back(i);
		}
	} else if ((k.y <= i.y) && (k.y <= j.y)) {
		points.push_back(k);
		if (i.y <= j.y) {
			points.push_back(i);
			points.push_back(j);		
		} else {
			points.push_back(j);
			points.push_back(i);
		}
	} 
	/*
	for (int i = 0; i < 3; i++) {
		std::cout << points[i].x << ", " << points[i].y << std::endl;
	}
	*/
	//intersect point
	//float intersect_y = points[1].y;
	float intersect_x = points[0].x + ((points[1].y - points[0].y) / (points[2].y - points[0].y)) * (points[2].x - points[0].x);
	//CanvasPoint intersect = CanvasPoint(intersect_x, intersect_y);
	//std::cout << "Intersect point: " << intersect.x << ", " << intersect.y << std::endl;
	//drawLine(window, points[1], intersect, c);	
	//fill top
	float topHeight = points[1].y - points[0].y;
	float boundary1StepSize = (intersect_x - points[0].x) / topHeight;
	float boundary2StepSize = (points[1].x - points[0].x) / topHeight;
	//std::cout << "Top triangle height: " << topHeight << std::endl;
	//std::cout << "Boundary 1 step size: " << boundary1StepSize << std::endl;
	//std::cout << "Boundary 2 step size: " << boundary2StepSize << std::endl;
	for (float row = points[0].y; row < points[1].y; row++) {
		//std::cout << "Row: " << row << std::endl;
		//boundaries
		float boundary1 = points[0].x + ((row - points[0].y) * boundary1StepSize);
		float boundary2 = points[0].x + ((row - points[0].y) * boundary2StepSize);
		//std::cout << "Row: " << row << std::endl << "Boundaries: " << boundary1 << ", " << boundary2;
		CanvasPoint p1 = CanvasPoint(std::round(boundary1), row);
		CanvasPoint p2 = CanvasPoint(std::round(boundary2), row);
		drawLine(window, p1, p2, c);
	}	
	//fill bottom
	float bottomHeight = points[2].y - points[1].y;
	float boundary3StepSize = (points[2].x - intersect_x) / bottomHeight;
	float boundary4StepSize = (points[2].x - points[1].x) / bottomHeight;
	for (float row = points[1].y; row < points[2].y; row++) {
		float boundary3 = intersect_x + ((row - points[1].y) * boundary3StepSize);
		float boundary4 = points[1].x + ((row - points[1].y) * boundary4StepSize);
		CanvasPoint p3 = CanvasPoint(std::round(boundary3), row);
		CanvasPoint p4 = CanvasPoint(std::round(boundary4), row);
		drawLine(window, p3, p4, c);
	}
	//draw white border
	Colour w = Colour(255,255,255);
	drawLine(window, i, j, w);
	drawLine(window, i, k, w);
	drawLine(window, j, k, w);
	
}

/*
int distance(CanvasPoint a, CanvasPoint b) {
	return 
}
*/

void mapTexture(DrawingWindow &window, CanvasTriangle t, TextureMap m) {
	//std::cout << "Height: " << mapWidth << std::endl;
	//std::cout << "Width: " << mapHeight << std::endl;
	
	std::vector<CanvasPoint> trianglePoints;
	std::vector<TexturePoint> texturePoints;
	CanvasPoint i = t.vertices[0];
	CanvasPoint j = t.vertices[1];
	CanvasPoint k = t.vertices[2];	
	if ((i.y <= j.y) && (i.y <= k.y)) {
		trianglePoints.push_back(i);
		texturePoints.push_back(i.texturePoint);
		if (j.y <= k.y) {
			trianglePoints.push_back(j);
			texturePoints.push_back(j.texturePoint);
			trianglePoints.push_back(k);
			texturePoints.push_back(k.texturePoint);
		} else {
			trianglePoints.push_back(k);
			texturePoints.push_back(k.texturePoint);
			trianglePoints.push_back(j);
			texturePoints.push_back(j.texturePoint);
		}
	} else if ((j.y <= i.y) && (j.y <= k.y)) {
		trianglePoints.push_back(j);
		texturePoints.push_back(j.texturePoint);
		if (i.y <= k.y) {
			trianglePoints.push_back(i);
			texturePoints.push_back(i.texturePoint);
			trianglePoints.push_back(k);
			texturePoints.push_back(k.texturePoint);
		} else {
			trianglePoints.push_back(k);
			texturePoints.push_back(k.texturePoint);
			trianglePoints.push_back(i);
			texturePoints.push_back(i.texturePoint);
		}
	} else if ((k.y <= i.y) && (k.y <= j.y)) {
		trianglePoints.push_back(k);
		texturePoints.push_back(k.texturePoint);
		if (i.y <= j.y) {
			trianglePoints.push_back(i);
			texturePoints.push_back(i.texturePoint);
			trianglePoints.push_back(j);		
			texturePoints.push_back(j.texturePoint);
		} else {
			trianglePoints.push_back(j);
			texturePoints.push_back(j.texturePoint);
			trianglePoints.push_back(i);
			texturePoints.push_back(i.texturePoint);
		}
	} 
	/*
	for (int i = 0; i < 3; i++) {
		std::cout << trianglePoints[i].x << ", " << trianglePoints[i].y << std::endl;
	}
	for (int i = 0; i < 3; i++) {
		std::cout << texturePoints[i].x << ", " << texturePoints[i].y << std::endl;
	}
	Colour w = Colour(255,255,255);
	drawLine(window, i, j, w);
	drawLine(window, i, k, w);
	drawLine(window, j, k, w);
	*/
	float intersectRatio = ((trianglePoints[1].y - trianglePoints[0].y) / (trianglePoints[2].y - trianglePoints[0].y));
	//std::cout << "intersect ratio: " << intersectRatio << std::endl;
	// intersect canvas point
	//(y1-y0) / (y2-y0) * (x2 - x0)
	float intersect_x = trianglePoints[0].x + intersectRatio * (trianglePoints[2].x - trianglePoints[0].x);
	
	
	//CanvasPoint intersect = CanvasPoint(intersect_x, trianglePoints[1].y);
	//drawLine(window, trianglePoints[1], intersect, w);
	
	// intersect texture point
	float intersect_texture_x = texturePoints[0].x + intersectRatio * (texturePoints[2].x - texturePoints[0].x);
	float intersect_texture_y = texturePoints[0].y + intersectRatio * (texturePoints[2].y - texturePoints[0].y);
	
	//std::cout << "Intersect point: " << intersect_texture_x << ", " << intersect_texture_y << std::endl;
		
	// fill top
	//number of rows in top triangle
	float rows = trianglePoints[1].y - trianglePoints[0].y;
	float b_p0p1_s = (trianglePoints[1].x - trianglePoints[0].x) / rows;
	float b_p0pi_s = (intersect_x - trianglePoints[0].x) / rows;
	float b_p0p1_m_x_s = (texturePoints[1].x - texturePoints[0].x) / rows;
	float b_p0p1_m_y_s = (texturePoints[1].y - texturePoints[0].y) / rows;
	float b_p0pi_m_x_s = (intersect_texture_x - texturePoints[0].x) / rows;
	float b_p0pi_m_y_s = (intersect_texture_y - texturePoints[0].y) / rows;
	for (float row = trianglePoints[0].y; row < trianglePoints[1].y; row++) {
		float b_p0p1_m_x = texturePoints[0].x + ((row - trianglePoints[0].y) * b_p0p1_m_x_s);
		float b_p0p1_m_y = texturePoints[0].y + ((row - trianglePoints[0].y) * b_p0p1_m_y_s);
		float b_p0pi_m_x = texturePoints[0].x + ((row - trianglePoints[0].y) * b_p0pi_m_x_s);
		float b_p0pi_m_y = texturePoints[0].y + ((row - trianglePoints[0].y) * b_p0pi_m_y_s);
		float b_p0p1 = trianglePoints[0].x + ((row - trianglePoints[0].y) * b_p0p1_s);
		float b_p0pi = trianglePoints[0].x + ((row - trianglePoints[0].y) * b_p0pi_s);
		float start = std::min(b_p0p1, b_p0pi);
		float end = std::max(b_p0p1, b_p0pi);
		float pointCount = std::abs(end - start);
		float step_size = (end - start) / pointCount;
		float start_m_x;
		float start_m_y;
		float end_m_x;
		float end_m_y;
		if (start == b_p0p1) {
			start_m_x = b_p0p1_m_x;
			start_m_y = b_p0p1_m_y;
			end_m_x = b_p0pi_m_x;
			end_m_y = b_p0pi_m_y;
		} else {
			start_m_x = b_p0pi_m_x;
			start_m_y = b_p0pi_m_y;
			end_m_x = b_p0p1_m_x;
			end_m_y = b_p0p1_m_y;
		}
		float step_size_m_x = (end_m_x - start_m_x) / pointCount;
		float step_size_m_y = (end_m_y - start_m_y) / pointCount;		
		for (float point = start; point < end; point++) {
			float m_x = start_m_x + ((point - start) * step_size_m_x);
			float m_y = start_m_y + ((point - start) * step_size_m_y);
			window.setPixelColour(std::round(point), std::round(row), m.pixels[(std::round(m_y) * m.width) + m_x]);
		}
	} 
	//fill bottom
	rows = trianglePoints[2].y - trianglePoints[1].y;
	float b_p1p2_s = (trianglePoints[2].x - trianglePoints[1].x) / rows;
	float b_pip2_s = (trianglePoints[2].x - intersect_x) / rows;
	float b_p1p2_m_x_s = (texturePoints[2].x - texturePoints[1].x) / rows;
	float b_p1p2_m_y_s = (texturePoints[2].y - texturePoints[1].y) / rows;
	float b_pip2_m_x_s = (texturePoints[2].x - intersect_texture_x) / rows;
	float b_pip2_m_y_s = (texturePoints[2].y - intersect_texture_y) / rows;
	for (float row = trianglePoints[1].y; row < trianglePoints[2].y; row++) {
		float b_p1p2_m_x = texturePoints[1].x + ((row - trianglePoints[1].y) * b_p1p2_m_x_s);
		float b_p1p2_m_y = texturePoints[1].y + ((row - trianglePoints[1].y) * b_p1p2_m_y_s);
		float b_pip2_m_x = intersect_texture_x + ((row - trianglePoints[1].y) * b_pip2_m_x_s);
		float b_pip2_m_y = intersect_texture_y + ((row - trianglePoints[1].y) * b_pip2_m_y_s);
		float b_p1p2 = trianglePoints[1].x + ((row - trianglePoints[1].y) * b_p1p2_s);
		float b_pip2 = intersect_x + ((row - trianglePoints[1].y) * b_pip2_s);
		float start = std::min(b_p1p2, b_pip2);
		float end = std::max(b_p1p2, b_pip2);
		float pointCount = std::abs(end - start);
		float step_size = (end - start) / pointCount;
		float start_m_x;
		float start_m_y;
		float end_m_x;
		float end_m_y;
		if (start == b_p1p2) {
			start_m_x = b_p1p2_m_x;
			start_m_y = b_p1p2_m_y;
			end_m_x = b_pip2_m_x;
			end_m_y = b_pip2_m_y;
		} else {
			start_m_x = b_pip2_m_x;
			start_m_y = b_pip2_m_y;
			end_m_x = b_p1p2_m_x;
			end_m_y = b_p1p2_m_y;
		}
		float step_size_m_x = (end_m_x - start_m_x) / pointCount;
		float step_size_m_y = (end_m_y - start_m_y) / pointCount;		
		for (float point = start; point < end; point++) {
			float m_x = start_m_x + ((point - start) * step_size_m_x);
			float m_y = start_m_y + ((point - start) * step_size_m_y);
			window.setPixelColour(std::round(point), std::round(row), m.pixels[(std::round(m_y) * m.width) + m_x]);
		}
	}
	
}

void map(DrawingWindow &window, TextureMap t) {
	for (size_t y = 0; y < t.height; y++) {
		for (size_t x = 0; x < t.width; x++) {
			window.setPixelColour(x, y, t.pixels[y*t.width+ x]);
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) drawStrokedTriangle(window, CanvasTriangle(CanvasPoint(rand() % window.width, rand() % window.height), CanvasPoint(rand() % window.width, rand() % window.height), CanvasPoint(rand() % window.width, rand() % window.height)), Colour(rand() % 256, rand() % 256, rand() % 256));
		else if (event.key.keysym.sym == SDLK_f) drawFilledTriangle(window, CanvasTriangle(CanvasPoint(rand() % window.width, rand() % window.height), CanvasPoint(rand() % window.width, rand() % window.height), CanvasPoint(rand() % window.width, rand() % window.height)), Colour(rand() % 256, rand() % 256, rand() % 256));
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	srand(time(NULL));
	// Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
	// Colour d = Colour(rand() % 256, rand() % 256, rand() % 256);
	// Colour e = Colour(rand() % 256, rand() % 256, rand() % 256);
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		CanvasPoint p1 = CanvasPoint(160, 10);
		p1.texturePoint = TexturePoint(195, 5);
		CanvasPoint p2 = CanvasPoint(300, 230);
		p2.texturePoint = TexturePoint(395,380);
		CanvasPoint p3 = CanvasPoint(10, 150);
		p3.texturePoint = TexturePoint(65, 330);
		mapTexture(window, CanvasTriangle(p1, p2, p3), TextureMap("texture.ppm"));
		//map(window, TextureMap("texture.ppm"));
		// Colour c = Colour(rand() % 256, rand() % 256, rand() % 256);
		// drawLine(window, CanvasPoint(0.0,0.0), CanvasPoint(window.width/2, window.height/2), c);
		// Colour d = Colour(rand() % 256, rand() % 256, rand() % 256);
		// drawLine(window, CanvasPoint(window.width/2, 0), CanvasPoint(window.width/2, window.height), d);
		// Colour e = Colour(rand() % 256, rand() % 256, rand() % 256);
		// drawLine(window, CanvasPoint(window.width/6, window.height/2), CanvasPoint(window.width*5/6, window.height/2), e);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
