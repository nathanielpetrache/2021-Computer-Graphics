#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <ModelTriangle.h>
#include <glm/glm.hpp>
#include <iomanip>
#include <Colour.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <CanvasPoint.h>
#include <limits>


#define WIDTH 1000
#define HEIGHT 1000
#define PI M_PI

glm::vec3 cameraPosition(0.0, 0.0, 4.0);
float focalLength = 2.0;
glm::mat3 cameraOrientation(
	glm::vec3(1.0, 0.0, 0.0),
	glm::vec3(0.0, 1.0, 0.0),
	glm::vec3(0.0, 0.0, 1.0)
);

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

std::vector<Colour> parseMTL(std::string filename) {
	std::vector<Colour> palette;
	//std::cout << "parsing MTL" << "\n";
	std::ifstream mtlFile(filename);
	if (mtlFile.is_open()) {
		//std::cout << "MTL file open.\n";
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
		//std::cout << "OBJ file open.\n";
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

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition) {
	//transpose from model coordinate system to camera coordinate system and scale
	float xyScale = HEIGHT/2;
	//vertexPosition = vertexPosition * cameraOrientation;
	glm::vec3 vertexPositionFromCamera(-(xyScale)*(vertexPosition[0]-cameraPosition[0]), (xyScale)*(vertexPosition[1]-cameraPosition[1]), (vertexPosition[2]-cameraPosition[2]));
	vertexPositionFromCamera = vertexPositionFromCamera * cameraOrientation;
	float canvasXFromCamera = (focalLength * (vertexPositionFromCamera[0] / vertexPositionFromCamera[2])) + (WIDTH / 2);
	float canvasYFromCamera = (focalLength * (vertexPositionFromCamera[1] / vertexPositionFromCamera[2])) + (HEIGHT / 2);
	return CanvasPoint(canvasXFromCamera, canvasYFromCamera, vertexPositionFromCamera[2]);
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

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c) {
	//window.clearPixels()
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

std::vector<CanvasPoint> sortPoints(CanvasPoint i, CanvasPoint j, CanvasPoint k){
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
	return points;
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle t, Colour c) {
	std::vector<CanvasPoint> points = sortPoints(t.vertices[0], t.vertices[1], t.vertices[2]);
	//intersect point
	//float intersect_y = points[1].y;
	float intersect_x = points[0].x + ((points[1].y - points[0].y) / (points[2].y - points[0].y)) * (points[2].x - points[0].x);
	//fill top
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
	//Colour w = Colour(255,255,255);
	
	drawLine(window, t.vertices[0], t.vertices[1], c);
	drawLine(window, t.vertices[0], t.vertices[2], c);
	drawLine(window, t.vertices[1], t.vertices[2], c);	
}

std::vector<std::vector<float>> drawLineWithDepth(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour c, std::vector<std::vector<float>> depthBuffer) {
	uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + (c.blue);
	//get x difference
	float xDiff = to.x - from.x;
	//get y difference
	float yDiff = to.y - from.y;
	//get z difference
	float zDiff = to.depth - from.depth;
	//get number of steps (max(x difference, y difference))
	float numberOfSteps = std::max(std::abs(xDiff), std::abs(yDiff));
	//get x step size
	float xStepSize = xDiff / numberOfSteps;
	//get y step size
	float yStepSize = yDiff / numberOfSteps;
	//get z step size
	float zStepSize = zDiff / numberOfSteps;
	//for loop
	//for each step in number of steps
	for (float step = 0.0; step < numberOfSteps; step++) {
		//get x
		float x = from.x + (step * xStepSize);
		//get y
		float y = from.y + (step * yStepSize);
		//get z
		float z = from.depth + (step * zStepSize);
		//z = ()
		//get 1/z of current point
		float currentPoint_1_Over_Z = 1 / z;
		//get 1/z value from depthBuffer
		float dB_1_Over_Z = depthBuffer[std::round(x)][std::round(y)];
		if (currentPoint_1_Over_Z < dB_1_Over_Z) {
			depthBuffer[std::round(x)][std::round(y)] = currentPoint_1_Over_Z;
			window.setPixelColour(std::round(x), std::round(y), colour);
		}
	}
	// return updated depth buffer
	return depthBuffer;
}



std::vector<std::vector<float>> drawFilledTriangleWithDepth(DrawingWindow &window, CanvasTriangle t, Colour c, std::vector<std::vector<float>> depthBuffer) {
	//sort i,j,k vertices by ascending y coordinates 
	std::vector<CanvasPoint> canvas_points = sortPoints(t.vertices[0], t.vertices[1], t.vertices[2]);
	std::vector<glm::vec3> points;
	for (int current_vec = 0; current_vec < 3; current_vec++) {
		points.push_back(glm::vec3(canvas_points[current_vec].x, canvas_points[current_vec].y, canvas_points[current_vec].depth));
	}
	//get intersect point: intersect x (forms line with j that horizontally splits triangle)
	float intersect_x = points[0][0] + ((points[1][1] - points[0][1]) / (points[2][1] - points[0][1])) * (points[2][0] - points[0][0]);
	//get intersect point depth (z from camera)
	float intersect_z = points[0][2] + ((points[1][1] - points[0][1]) / (points[2][1] - points[0][1])) * (points[2][2] - points[0][2]);

	glm::vec3 intersect_point = glm::vec3(intersect_x, points[1][1], intersect_z);
	//fill top
	//calculate top triangle height
	float topHeight = points[1][1] - points[0][1];

	std::vector<glm::vec3> left_boundary_line = interpolateThreeElementValues(points[0], intersect_point, topHeight+2);
	std::vector<glm::vec3> right_boundary_line = interpolateThreeElementValues(points[0], points[1], topHeight+2);

	for (int row_index = 0; row_index < topHeight+1; row_index++) {
		int row_length = std::abs(std::round(left_boundary_line[row_index][0] - right_boundary_line[row_index][0]));
		
		//float row_length = left_boundary_line[row][0] - right_boundary_line[row][0];
		std::vector<glm::vec3>line=interpolateThreeElementValues(left_boundary_line[row_index],right_boundary_line[row_index],row_length+2);
		for (int point_index=0;point_index<row_length+2;point_index++){
			float point_1_over_z=1.0/line[point_index][2];
			if(point_1_over_z<depthBuffer[std::round(line[point_index][0])][std::round(line[point_index][1])]){
				depthBuffer[std::round(line[point_index][0])][std::round(line[point_index][1])]=point_1_over_z;
				uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + (c.blue);
				window.setPixelColour(std::round(line[point_index][0]),std::round(line[point_index][1]),colour);
			}
		}

	}

	//fill bottom
	float bottomHeight=points[2][1]-points[1][1];
	left_boundary_line = interpolateThreeElementValues(intersect_point,points[2],bottomHeight+2);
	right_boundary_line = interpolateThreeElementValues(points[1],points[2],bottomHeight+2);
	for (int row_index = 0; row_index < bottomHeight+1; row_index++) {
		int row_length = std::abs(std::round(left_boundary_line[row_index][0] - right_boundary_line[row_index][0]));
		
		//float row_length = left_boundary_line[row][0] - right_boundary_line[row][0];
		std::vector<glm::vec3>line=interpolateThreeElementValues(left_boundary_line[row_index],right_boundary_line[row_index],row_length+2);
		for (int point_index=0;point_index<row_length+2;point_index++){
			float point_1_over_z=1.0/line[point_index][2];
			if(point_1_over_z<depthBuffer[std::round(line[point_index][0])][std::round(line[point_index][1])]){
				depthBuffer[std::round(line[point_index][0])][std::round(line[point_index][1])]=point_1_over_z;
				uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + (c.blue);
				window.setPixelColour(std::round(line[point_index][0]),std::round(line[point_index][1]),colour);
			}
		}

	}

	return depthBuffer;

	/*
	//calculate boundary step size for i.x to intersect x and i.x to j.x  
	float boundary1StepSize = (intersect_x - points[0].x) / topHeight;
	float boundary2StepSize = (points[1].x - points[0].x) / topHeight;
	//calculate z step size from i to intersect and i to j
	float boundary1_Z_StepSize = (intersect_z - points[0].depth) / topHeight;
	float boundary2_Z_StepSize = (points[1].depth - points[0].depth) / topHeight;
	//for loop
	//for each y value between i.y to j.y/intersect.y
	for (float row = points[0].y; row < points[1].y; row++) {
		//boundaries
		//get boundary x values
		float boundary1 = points[0].x + ((row - points[0].y) * boundary1StepSize);
		float boundary2 = points[0].x + ((row - points[0].y) * boundary2StepSize);
		//get boundary z values
		float boundary1_z = points[0].depth + ((row - points[0].y) * boundary1_Z_StepSize);
		float boundary2_z = points[0].depth + ((row - points[0].y) * boundary2_Z_StepSize);
		//get the two points of form: (boundary x, current y, z)
		CanvasPoint p1 = CanvasPoint(std::round(boundary1), row, boundary1_z);
		CanvasPoint p2 = CanvasPoint(std::round(boundary2), row, boundary2_z);
		//draw line between two points, passing in depth buffer, passing in model triangle
		depthBuffer = drawLineWithDepth(window, p1, p2, c, depthBuffer);
	}	
	//fill bottom
	//draw line between two points, passing in depth buffer, passing in model triangle
	//calculate bottom triangle height
	float bottomHeight = points[2].y - points[1].y;
	//calculate boundary step size for intersect x to k.x and j.x to k.x
	float boundary3StepSize = (points[2].x - intersect_x) / bottomHeight;
	float boundary4StepSize = (points[2].x - points[1].x) / bottomHeight;
	//calculate z step size from intersect to k and from j to k
	float boundary3_Z_StepSize = (points[2].depth - intersect_z) / bottomHeight;
	float boundary4_Z_StepSize = (points[2].depth - points[1].depth) / bottomHeight;
	//for loop
	//for each y value between j.y to k.y
	for (float row = points[1].y; row < points[2].y; row++) {
		//get boundary x values
		float boundary3 = intersect_x + ((row - points[1].y) * boundary3StepSize);
		float boundary4 = points[1].x + ((row - points[1].y) * boundary4StepSize);
		//get boundary z values
		float boundary3_z = intersect_z + ((row - points[1].y) * boundary3_Z_StepSize);
		float boundary4_z = intersect_z + ((row - points[1].y) * boundary4_Z_StepSize);
		//get the two points of form: (boundary x, current y, z)
		CanvasPoint p3 = CanvasPoint(std::round(boundary3), row, boundary3_z);
		CanvasPoint p4 = CanvasPoint(std::round(boundary4), row, boundary4_z);
		depthBuffer = drawLineWithDepth(window, p3, p4, c, depthBuffer);
	}
	
	depthBuffer = drawLineWithDepth(window, t.vertices[0], t.vertices[1], c, depthBuffer);
	depthBuffer = drawLineWithDepth(window, t.vertices[0], t.vertices[2], c, depthBuffer);
	depthBuffer = drawLineWithDepth(window, t.vertices[1], t.vertices[2], c, depthBuffer);
	
	return depthBuffer;
	*/
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

std::vector<std::vector<float>> initialiseDepthBuffer(int width, int height){
	std::vector<std::vector<float>> depthBuffer;
	//initialise depth buffer
	for (int x = 0; x < WIDTH; x++) {
		std::vector<float> row;
		for (int y = 0; y < HEIGHT; y++) {
			row.push_back(1.0/std::numeric_limits<float>::max());
		}
		depthBuffer.push_back(row);
	}
	return depthBuffer;
}

void drawRasterisedRenderWithDepth(DrawingWindow &window, std::vector<ModelTriangle> modelTriangles, std::vector<CanvasTriangle> canvasTriangles) {
	std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
	//for each canvas triangle
	for (int cT = 0; cT < canvasTriangles.size(); cT++) {
	//draw filled triangle (pass depth buffer)
		depthBuffer = drawFilledTriangleWithDepth(window, canvasTriangles[cT], modelTriangles[cT].colour, depthBuffer);
	}
}

void draw(DrawingWindow &window) {
	window.clearPixels();
	std::vector<ModelTriangle> triangles = parseOBJ("cornell-box.obj");
	std::vector<CanvasPoint> projectedPoints;
	std::vector<CanvasTriangle> canvasTriangles;
	std::vector<Colour> colourVector;
	//float focalLength = 2.0;
	for (int t = 0; t < triangles.size(); t++) {
		CanvasTriangle cT = CanvasTriangle();
		for (int v = 0; v < 3; v++) {
			//std::cout << "focal length = " << focalLength << std::endl;
			CanvasPoint projectedPoint = getCanvasIntersectionPoint(triangles[t].vertices[v]);
			//focalLength = cameraPosition[2]/2;
			//std::cout << "focal length 2 = " << focalLength << std::endl;
			projectedPoints.push_back(projectedPoint);
			cT.vertices[v] = projectedPoint;
		}
		canvasTriangles.push_back(cT);
	}
	
	//drawPointcloud(window, projectedPoints);
	//drawWireframe(window, canvasTriangles);
	//drawRasterisedRender(window, triangles, canvasTriangles);
	drawRasterisedRenderWithDepth(window, triangles, canvasTriangles);	
	
}

glm::mat3 xRotation(float theta){
	glm::mat3 mat = glm::mat3(
		glm::vec3(1.0, 0.0, 0.0), 
		glm::vec3(0.0, cos(theta), sin(theta)), 
		glm::vec3(0.0, -sin(theta), cos(theta))
	);
	return mat;
}

void rotateAboutX(float theta) {
	glm::mat3 mat = xRotation(theta);
	cameraPosition = mat * cameraPosition;
	cameraOrientation = mat * cameraOrientation;
}

glm::mat3 yRotation(float theta){
	glm::mat3 mat = glm::mat3(
		glm::vec3(cos(theta), 0.0, -sin(theta)), 
		glm::vec3(0.0, 1.0, 0.0), 
		glm::vec3(sin(theta), 0.0, cos(theta))
	);
	return mat;
}

void rotateAboutY(float theta) {
	glm::mat3 mat = yRotation(theta);
	cameraPosition = mat * cameraPosition;
	//std::cout << cameraPosition[0] << ", " << cameraPosition[1] << ", " << cameraPosition[2] << std::endl;
	//std::cout << "Camera Orientation: " << std::endl;
	//std::cout << cameraOrientation[0][0] << ", " << cameraOrientation[0][1] << ", " << cameraOrientation[0][2] << std::endl;
	//std::cout << cameraOrientation[1][0] << ", " << cameraOrientation[1][1] << ", " << cameraOrientation[1][2] << std::endl;
	//std::cout << cameraOrientation[2][0] << ", " << cameraOrientation[2][1] << ", " << cameraOrientation[2][2] << std::endl;
	cameraOrientation = mat * cameraOrientation;
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		
		if (event.key.keysym.sym == SDLK_LEFT) {
			std::cout << "LEFT" << std::endl;
			cameraPosition[0] += 0.05;
		} else if (event.key.keysym.sym == SDLK_RIGHT) {
			std::cout << "RIGHT" << std::endl;
			cameraPosition[0] -= 0.05;
		} else if (event.key.keysym.sym == SDLK_UP) {
			std::cout << "UP" << std::endl;
			cameraPosition[1] -= 0.05;
		} else if (event.key.keysym.sym == SDLK_DOWN) {
			std::cout << "DOWN" << std::endl;
			cameraPosition[1] += 0.05;
		} else if (event.key.keysym.sym == SDLK_PAGEUP) {
			std::cout << "ZOOM OUT" << std::endl;
			cameraPosition[2] += 0.1;
		} else if (event.key.keysym.sym == SDLK_PAGEDOWN) {
			std::cout << "ZOOM IN" << std::endl;
			cameraPosition[2] -= 0.1;			
		} else if (event.key.keysym.sym == SDLK_w) {
			std::cout << "X ROTATE" << std::endl;
			rotateAboutX(0.01);
		} else if (event.key.keysym.sym == SDLK_s) {
			std::cout << "X ROTATE" << std::endl;
			rotateAboutX(-0.01);
		}else if (event.key.keysym.sym == SDLK_a) {
			std::cout << "Y ROTATE" << std::endl;
			rotateAboutY(0.01);
		} else if (event.key.keysym.sym == SDLK_d) {
			std::cout << "Y ROTATE" << std::endl;
			rotateAboutY(-0.01);
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	//return cameraPosition;
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) {
			//updates camera position + orientation
			handleEvent(event, window);		
		}
		window.clearPixels();
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
