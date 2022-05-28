#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <ModelTriangle.h>
#include <glm/glm.hpp>
#include <Colour.h>
#include <CanvasPoint.h>
#include <cmath>

#define WIDTH 320
#define HEIGHT 240

glm::vec3 cameraPosition(0.0, 0.0, 4.0);

std::vector<Colour> parseMTL(std::string filename) {
	std::vector<Colour> palette;
	std::cout << "parsing MTL" << "\n";
	std::ifstream mtlFile(filename);
	if (mtlFile.is_open()) {
		std::cout << "MTL file open.\n";
		std::string line;
		std::string currentColourName;
		//each colour stored as 0-255 in Colour object
		std::vector<int> rgb;
		//colour count for printing colour from vector
		int colourIndex = 0;
		while (getline(mtlFile, line)){
			std::string delimiter = " ";
			//get start token
			std::string token = line.substr(0, line.find(delimiter));
			//remove token from line
			line.erase(0, line.find(delimiter)+1);
			if (token == "newmtl") {
				//colour name is remainder of line
				currentColourName = line;
				//get next line
				getline(mtlFile, line);
				std::string kd_token = "Kd";
				//remove kd token
				line.erase(0, line.find(kd_token)+kd_token.length()+1);
				float currentColourComponent;
				for (int i = 0; i < 3; i++) {
					//get current float value from line
					currentColourComponent = std::stof(line.substr(0, line.find(delimiter)-1));

					rgb.push_back(round(255*currentColourComponent));
					line.erase(0, line.find(delimiter)+1);
				}			
				Colour newColour = Colour(currentColourName, rgb[3*colourIndex], rgb[3*colourIndex+1], rgb[3*colourIndex+2]);
				palette.push_back(newColour);
				colourIndex++;
			}
		}
		mtlFile.close();
 	} else {
		std::cout << "MTL file not open.\n";
	}
	return palette;
}

std::vector<ModelTriangle> parseOBJ(std::string filename){
	std::vector<glm::vec3> points;
	std::vector<ModelTriangle> faces;
	std::vector<Colour> palette;
	float scalingFactor = 0.17;
	//file as input stream
	std::ifstream objFile(filename);
	//check if open
	if (objFile.is_open()) {
		std::cout << "OBJ file open.\n";
		std::string line;
		std::string colourName;
		while (getline(objFile, line)){
			//parse line, delimeter = " "
			std::string delimiter = " ";
			//get start token
			std::string token = line.substr(0, line.find(delimiter));
			//remove token from line
			line.erase(0, line.find(delimiter)+1);
			//depending on token, parse line
			if (token == "mtllib") {
				//get filename
				std::string mtlFile = line;
				//parse file
				palette = parseMTL(mtlFile);
			} else if (token == "usemtl") {
				//get colour, update current colour name
				colourName = line;
			} else if (token == "v") {
				delimiter = " ";
				std::vector<float> currentPoint;
				float currentCoordinate;
				for (int i = 0; i < 3; i++) {
					currentCoordinate = std::stof(line.substr(0, line.find(delimiter)-1));
					currentCoordinate *= scalingFactor;
					currentPoint.push_back(currentCoordinate);
					line.erase(0, line.find(delimiter)+1);
				}
				glm::vec3 currentPointVec3(currentPoint[0], currentPoint[1], currentPoint[2]);
				points.push_back(currentPointVec3);
			} else if (token == "f") {
				delimiter = "/";
				int vertexIndex;
				ModelTriangle currentTriangle = ModelTriangle();
				for (int i = 0; i < 3; i++) {
					vertexIndex = std::stoi(line.substr(0, line.find(delimiter)))-1;
					currentTriangle.vertices[i] = points[vertexIndex];
					line.erase(0, line.find(delimiter)+1);
				}
				//apply current colour to triangle
				int j = 0;
				while (palette[j].name != colourName && j < palette.size()) {
					j++;
				}
				currentTriangle.colour = palette[j];
				faces.push_back(currentTriangle);
			}
		}
		objFile.close();
	} else {
		std::cout << "OBJ file not open.\n";
	}
	return faces;
}

std::vector<std::vector<float>> initializeDepthBuffer(int width, int height){
	std::vector<std::vector<float>> depthBuffer;
	for (int x = 0; x < width; x++) {
		std::vector<float> row;
		for (int y = 0; y < height; y++) {
			row.push_back(std::numeric_limits<float>::max());
		}
		depthBuffer.push_back(row);
	}
	return depthBuffer;
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

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float focalLength) {
	//transpose from model coordinate system to camera coordinate system and scale
	float xyScale = 250;
	glm::vec3 vertexPositionFromCamera(-(xyScale)*(vertexPosition[0]-cameraPosition[0]), (xyScale)*(vertexPosition[1]-cameraPosition[1]), vertexPosition[2]-cameraPosition[2]);
	//std::cout << "vertex position = (" << vertexPosition[0] << ", " << vertexPosition[1] << ", " << vertexPosition[2] << ")" << "\n";
	//std::cout << "vertex position from camera = (" << vertexPositionFromCamera[0] << ", " << vertexPositionFromCamera[1] << ", " << vertexPositionFromCamera[2] << ")" << "\n";
	float canvasXFromCamera = (focalLength * (vertexPositionFromCamera[0] / vertexPositionFromCamera[2])) + (WIDTH / 2);
	float canvasYFromCamera = (focalLength * (vertexPositionFromCamera[1] / vertexPositionFromCamera[2])) + (HEIGHT / 2);
	float zDepth = sqrt(pow(vertexPositionFromCamera[0],2) + pow(vertexPositionFromCamera[1], 2) + pow(vertexPositionFromCamera[2], 2));
	return CanvasPoint(canvasXFromCamera, canvasYFromCamera, zDepth);
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

std::vector<CanvasPoint> sortTrianglePoints(CanvasPoint i, CanvasPoint j, CanvasPoint k) {
	std::vector<CanvasPoint> points;
	if ((i.y <= j.y) && (i.y <= k.y)) {
		points.push_back(i);
		if (j.y <= k.y) {
			points.push_back(j);
			points.push_back(k);
		} 
		else {
			points.push_back(k);
			points.push_back(j);
		}
	}
	else if ((j.y <= i.y) && (j.y <= k.y)) {
		points.push_back(j);
		if (i.y <= k.y) {
			points.push_back(i);
			points.push_back(k);
		} 
		else {
			points.push_back(k);
			points.push_back(i);
		}
	}
	else if ((k.y <= i.y) && (k.y <= j.y)) {
		points.push_back(k);
		if (i.y <= j.y) {
			points.push_back(i);
			points.push_back(j);		
		}
		else {
			points.push_back(j);
			points.push_back(i);
		}
	} 
	return points;
}

void fillTopOfTriangle(DrawingWindow &window, std::vector<CanvasPoint> points, Colour c) {
	float intersect_x = points[0].x + ((points[1].y - points[0].y) / (points[2].y - points[0].y)) * (points[2].x - points[0].x);
	float topHeight = points[1].y - points[0].y;
	float boundary1StepSize = (intersect_x - points[0].x) / topHeight;
	float boundary2StepSize = (points[1].x - points[0].x) / topHeight;
	for (float row = points[0].y; row < points[1].y; row++) {
		//boundaries
		float boundary1 = points[0].x + ((row - points[0].y) * boundary1StepSize);
		float boundary2 = points[0].x + ((row - points[0].y) * boundary2StepSize);
		CanvasPoint p1 = CanvasPoint(std::round(boundary1), row);
		CanvasPoint p2 = CanvasPoint(std::round(boundary2), row);
		drawLine(window, p1, p2, c);
	}	
}

void fillBottomOfTriangle(DrawingWindow &window, std::vector<CanvasPoint> points, Colour c) {
	float intersect_x = points[0].x + ((points[1].y - points[0].y) / (points[2].y - points[0].y)) * (points[2].x - points[0].x);
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
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle t, Colour c) {
	std::vector<CanvasPoint> points = sortTrianglePoints(t.vertices[0], t.vertices[1], t.vertices[2]);
	fillTopOfTriangle(window, points, c);
	fillBottomOfTriangle(window, points, c);	
	drawLine(window, t.vertices[0], t.vertices[1], c);
	drawLine(window, t.vertices[0], t.vertices[2], c);
	drawLine(window, t.vertices[1], t.vertices[2], c);	
}

void drawPointcloud(DrawingWindow &window, std::vector<CanvasPoint> projectedPoints) {
	window.clearPixels();
	for (int i = 0; i < projectedPoints.size(); i++) {
		int red = 255;
		int green = 255;
		int blue = 255;
		uint32_t colour = (255 << 24) + (red << 16) + (green << 8) + blue;
		//projectedPoints[i].x *= 1.2;
		//projectedPoints[i].y *= 1.2;
		window.setPixelColour(projectedPoints[i].x, projectedPoints[i].y, colour);
	}
}

void drawWireframe(DrawingWindow &window, std::vector<CanvasTriangle> canvasTriangles) {
	Colour c = Colour(255,255,255);
	for (int t = 0; t<canvasTriangles.size(); t++) {
		drawStrokedTriangle(window, canvasTriangles[t], c);
	}
}

void drawRasterisedRender(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, std::vector<CanvasTriangle> canvasTriangles) {
	for (int cT = 0; cT < canvasTriangles.size(); cT++) {
		drawFilledTriangle(window, canvasTriangles[cT], modelTriangles[cT].colour);
	}
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	std::vector<ModelTriangle> triangles = parseOBJ("cornell-box.obj");
	std::vector<CanvasPoint> projectedPoints;
	std::vector<CanvasTriangle> canvasTriangles;
	std::vector<Colour> colourVector;
	float focalLength = 2.0;
	std::vector<std::vector<float>> depthBuffer = initializeDepthBuffer(WIDTH, HEIGHT);
	for (int t = 0; t < triangles.size(); t++) {
		CanvasTriangle cT = CanvasTriangle();
		for (int v = 0; v < 3; v++) {
			CanvasPoint projectedPoint = getCanvasIntersectionPoint(triangles[t].vertices[v], focalLength);
			projectedPoints.push_back(projectedPoint);
			cT.vertices[v] = projectedPoint;
		}
		canvasTriangles.push_back(cT);
	
	}
	//drawPointcloud(window, projectedPoints);
	//drawWireframe(window, canvasTriangles);
	drawRasterisedRender(window, triangles, canvasTriangles);
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
