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
std::vector<std::vector<float>> depthBuffer;
bool orbit_toggled_on = false;

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
	glm::vec3 vertexPositionFromCamera((vertexPosition[0]-cameraPosition[0]), (vertexPosition[1]-cameraPosition[1]), (vertexPosition[2]-cameraPosition[2]));
	vertexPositionFromCamera = vertexPositionFromCamera * cameraOrientation;
	float canvasXFromCamera = -xyScale*(focalLength * (vertexPositionFromCamera[0] / vertexPositionFromCamera[2])) + (WIDTH / 2);
	float canvasYFromCamera = xyScale*(focalLength * (vertexPositionFromCamera[1] / vertexPositionFromCamera[2])) + (HEIGHT / 2);
	return CanvasPoint(canvasXFromCamera, canvasYFromCamera, vertexPositionFromCamera[2]);
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
	
	drawLine(window, i, j, c);
	drawLine(window, i, k, c);
	drawLine(window, j, k, c);	
}

std::vector<std::vector<float>> fillHalfTriangle(DrawingWindow &window, glm::vec3 referencePoint, glm::vec3 p2, glm::vec3 p3, Colour c, std::vector<std::vector<float>> depthBuffer) {
	
	float height = std::abs(referencePoint[1] - p2[1]);

	std::vector<glm::vec3> left_boundary_line = interpolateThreeElementValues(referencePoint, p2, height+2);
	//std::cout<<"left boundary successful\n";
	std::vector<glm::vec3> right_boundary_line = interpolateThreeElementValues(referencePoint, p3, height+2);
	//std::cout<<"right boundary successful\n";
	for (int row_index = 0; row_index < height+1; row_index++) {
		int row_length = std::abs(std::round(left_boundary_line[row_index][0] - right_boundary_line[row_index][0]));
		//std::cout<<"row "<<row_index<<" length = "<<row_length<<"\n";
		std::vector<glm::vec3>line=interpolateThreeElementValues(left_boundary_line[row_index],right_boundary_line[row_index],row_length+2);
		//std::cout<<"line between boundaries successful\n";

		for (int point_index=0;point_index<row_length+2;point_index++){
			int x = std::round(line[point_index][0]);
			int y = std::round(line[point_index][1]);
			
			//std::cout<<"attempting to draw point ("<<x<<", "<<y<<")\n";
	
			if (((x>=0)&&(x<WIDTH))&&((y>=0)&&(y<HEIGHT))){
				//std::cout<<"point in bounds so continuing to draw point\n";
				float point_1_over_z=1.0/line[point_index][2];
				//std::cout<<"1/z value of point = "<<point_1_over_z<<"\n";
				//std::cout<<"1/z value of depthBuffer entry = "<<depthBuffer[x][y]<<"\n";
				if(point_1_over_z<depthBuffer[x][y]){
					//std::cout<<"1/z is less than depthBuffer entry so continuing to draw point\n";
					depthBuffer[x][y]=point_1_over_z;
					uint32_t colour = (255 << 24) + (c.red << 16) + (c.green << 8) + (c.blue);
					window.setPixelColour(x,y,colour);
				}
			}else{
				//std::cout<<"out of bounds\n";
			}
		}
	}
	return depthBuffer;
}
std::vector<std::vector<float>> drawFilledTriangleWithDepth(DrawingWindow &window, CanvasTriangle t, Colour c, std::vector<std::vector<float>> depthBuffer, int triangle_index) {
	//std::cout<<"drawing canvas triangle "<<triangle_index<<"\n";
	//std::cout<<"vertices = \n";
	//std::cout<<t.vertices[0]<<"\n";
	//std::cout<<t.vertices[1]<<"\n";
	//std::cout<<t.vertices[2]<<"\n";
	//sort i,j,k vertices by ascending y coordinates 
	std::vector<CanvasPoint> canvas_points = sortPoints(t.vertices[0], t.vertices[1], t.vertices[2]);
	std::vector<glm::vec3> points;
	for (int current_vec = 0; current_vec < 3; current_vec++) {
		points.push_back(glm::vec3(canvas_points[current_vec].x, canvas_points[current_vec].y, canvas_points[current_vec].depth));
	}
	float intersect_x = points[0][0] + ((points[1][1] - points[0][1]) / (points[2][1] - points[0][1])) * (points[2][0] - points[0][0]);
	float intersect_z = points[0][2] + ((points[1][1] - points[0][1]) / (points[2][1] - points[0][1])) * (points[2][2] - points[0][2]);
	glm::vec3 intersect_point = glm::vec3(intersect_x, points[1][1], intersect_z);
	// std::cout<<"drawing top half...\n";
	// std::cout<<"points = \n";
	// std::cout<<"("<<points[0][0]<<", "<<points[0][1]<<", "<<points[0][2]<<")\n";
	// std::cout<<"("<<points[1][0]<<", "<<points[1][1]<<", "<<points[1][2]<<")\n";
	// std::cout<<"("<<intersect_point[0]<<", "<<intersect_point[1]<<", "<<intersect_point[2]<<")\n";
	depthBuffer = fillHalfTriangle(window, points[0], points[1], intersect_point, c, depthBuffer);
	// std::cout<<"top half successfully drawn\n";
	// std::cout<<"drawing bottom half...\n";
	depthBuffer = fillHalfTriangle(window, points[2], points[1], intersect_point, c, depthBuffer);
	// std::cout<<"bottom half successfully drawn\n";
	return depthBuffer;
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

void drawRasterisedRenderWithDepth(DrawingWindow &window, std::vector<Colour> colours, std::vector<CanvasTriangle> canvasTriangles) {
	std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
	//for each canvas triangle
	for (int cT = 0; cT < canvasTriangles.size(); cT++) {
	//draw filled triangle (pass depth buffer)
		depthBuffer = drawFilledTriangleWithDepth(window, canvasTriangles[cT], colours[cT], depthBuffer, cT);
	}
}
void printDetails(float theta){
	std::cout<<"orientation angle = "<<theta<<"\n";
	std::cout << "camera position("<<cameraPosition[0] << ", " << cameraPosition[1] << ", " << cameraPosition[2]<<")\n";
	std::cout <<"camera orientation(\n";
	std::cout<<cameraOrientation[0][0]<<", "<<cameraOrientation[0][1]<<", "<<cameraOrientation[0][2]<<"\n";
	std::cout<<cameraOrientation[1][0]<<", "<<cameraOrientation[1][1]<<", "<<cameraOrientation[1][2]<<"\n";
	std::cout<<cameraOrientation[2][0]<<", "<<cameraOrientation[2][1]<<", "<<cameraOrientation[2][2]<<")\n";

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
	cameraPosition = xRotation(theta) * cameraPosition;
	//cameraOrientation = mat * cameraOrientation;
	printDetails(theta);
}
void orientAboutX(float theta){
	cameraOrientation = xRotation(theta) * cameraOrientation;
	printDetails(theta);
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
	cameraPosition = yRotation(theta) * cameraPosition;
	printDetails(theta);
}
void orientAboutY(float theta){
	cameraOrientation = yRotation(theta) * cameraOrientation;
	printDetails(theta);
}
void lookAt(glm::vec3 point){
	glm::vec3 forward = glm::normalize(glm::vec3(cameraPosition[0]-point[0],cameraPosition[1]-point[1],cameraPosition[2]-point[2]));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0.0, 1.0, 0.0),forward));
	glm::vec3 up = glm::cross(forward, right);

	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
}
void orbit(){
	float theta = PI /180;
	if (orbit_toggled_on) {
		cameraPosition = yRotation(theta) * cameraPosition;
		//cameraOrientation = yRotation(theta) * cameraOrientation;
		lookAt(glm::vec3(0.0,0.0,0.0));
		printDetails(theta);
	}
}

void draw(DrawingWindow &window, std::vector<ModelTriangle> triangles) {
	window.clearPixels();
	//std::vector<ModelTriangle> triangles = parseOBJ("cornell-box.obj");
	std::vector<CanvasPoint> projectedPoints;
	std::vector<CanvasTriangle> canvasTriangles;
	std::vector<Colour> colourVector;
	//float focalLength = 2.0;
	//rotateAboutY(PI/2);
	orbit();
	for (int t = 0; t < triangles.size(); t++) {
		CanvasTriangle cT = CanvasTriangle();
		int numberOfPointsOffCanvas = 0;
		//std::vector<CanvasPoint> currentTrianglePoints;
		for (int v = 0; v < 3; v++) {
			//std::cout << "focal length = " << focalLength << std::endl;
			CanvasPoint projectedPoint = getCanvasIntersectionPoint(triangles[t].vertices[v]);
			//focalLength = cameraPosition[2]/2;
			//std::cout << "focal length 2 = " << focalLength << std::endl;
			
			if ((projectedPoint.x >= WIDTH || projectedPoint.x < 0 || projectedPoint.y >= HEIGHT || projectedPoint.y < 0) || projectedPoint.depth > 0){
				//std::cout<<"point ("<<projectedPoint.x<<", "<<projectedPoint.y<<") is off canvas\n";
				numberOfPointsOffCanvas++;
				//currentTrianglePoints.push_back(projectedPoint);
			}
			cT.vertices[v] = projectedPoint;
		}
		if (numberOfPointsOffCanvas != 3){
			std::cout<<"triangle "<<cT<<" will be drawn\n";
			std::cout<<"("<<cT.vertices[0].x<<", "<<cT.vertices[0].y<<")\n";
			std::cout<<"("<<cT.vertices[1].x<<", "<<cT.vertices[1].y<<")\n";
			std::cout<<"("<<cT.vertices[2].x<<", "<<cT.vertices[2].y<<")\n";
			canvasTriangles.push_back(cT);
			colourVector.push_back(triangles[t].colour);
		}
	}
		//if (validPoints == 3) {
			//for (int point_index=0;point_index<3;point_index++){
				//CanvasPoint currentPoint = currentTrianglePoints[point_index];
				//projectedPoints.push_back(currentPoint);
			//}
		//canvasTriangles.push_back(cT);
		//colourVector.push_back(triangles[t].colour);
		//}
	std::cout<<"successful canvas triangles\n";
	std::cout<<"drawing "<<canvasTriangles.size()<<"/"<<triangles.size()<<" triangles\n";
	//drawWireframe(window, canvasTriangles);
	drawRasterisedRenderWithDepth(window, colourVector, canvasTriangles);	
}
void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {

		float theta = PI / 180;
		
		if (event.key.keysym.sym == SDLK_LEFT) {
			std::cout << "LEFT" << std::endl;
			cameraPosition[0] += 0.05;
		}else if (event.key.keysym.sym == SDLK_RIGHT) {
			std::cout << "RIGHT" << std::endl;
			cameraPosition[0] -= 0.05;
		}else if (event.key.keysym.sym == SDLK_UP) {
			std::cout << "UP" << std::endl;
			cameraPosition[1] -= 0.05;
		}else if (event.key.keysym.sym == SDLK_DOWN) {
			std::cout << "DOWN" << std::endl;
			cameraPosition[1] += 0.05;
		}else if (event.key.keysym.sym == SDLK_PAGEUP) {
			std::cout << "ZOOM OUT" << std::endl;
			cameraPosition[2] += 0.1;
		}else if (event.key.keysym.sym == SDLK_PAGEDOWN) {
			std::cout << "ZOOM IN" << std::endl;
			cameraPosition[2] -= 0.1;			
		}else if (event.key.keysym.sym == SDLK_w) {
			std::cout << "X ROTATE" << std::endl;
			rotateAboutX(-theta);
		}else if (event.key.keysym.sym == SDLK_s) {
			std::cout << "X ROTATE" << std::endl;
			rotateAboutX(theta);
		}else if (event.key.keysym.sym == SDLK_a) {
			std::cout << "Y ROTATE" << std::endl;
			rotateAboutY(-theta);
		}else if (event.key.keysym.sym == SDLK_d) {
			std::cout << "Y ROTATE" << std::endl;
			rotateAboutY(theta);
		}else if (event.key.keysym.sym == SDLK_i) {
			std::cout << "X ORIENT" << std::endl;
			orientAboutX(-theta);
		}else if (event.key.keysym.sym == SDLK_k) {
			std::cout << "X ORIENT" << std::endl;
			orientAboutX(theta);
		}else if (event.key.keysym.sym == SDLK_j) {
			std::cout << "Y ORIENT" << std::endl;
			orientAboutY(-theta);
		}else if (event.key.keysym.sym == SDLK_l) {
			std::cout << "Y ORIENT" << std::endl;
			orientAboutY(theta);
		}else if (event.key.keysym.sym == SDLK_o) {
			std::cout << "TOGGLED ORBIT" << std::endl;
			orbit_toggled_on = !orbit_toggled_on;
		}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
	//return cameraPosition;
}
int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	std::vector<ModelTriangle> triangles = parseOBJ("cornell-box.obj");
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) {
			//updates camera position + orientation
			handleEvent(event, window);		
		}
		//orientAboutY(PI/180);
		window.clearPixels();
		draw(window, triangles);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
