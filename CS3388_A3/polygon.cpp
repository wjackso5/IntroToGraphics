/*
Author: Will Jackson
Purpose: This is the implementation file for the class file that stores a polygon
Date: 11/8/2017
*/
#include "polygon.h"
//constructor for a polygon
//gets the list of vertices and then calculates the normals
Polygon::Polygon(list<dmatrix_t> vertexList){
	vertexList=vertexList;
	surfaceNormal = *dcross_product(*mat_sub(*vertexList[1], *vertexList[0]) , *mat_sub(*vertexList[2], *vertexList[1]))
}
//deconstuctor
Polygon::~Polygon(){
}
//getter for surface normal
dmatrix Polygon::getNormal(){
	return surfaceNormal;
}
void Polygon::fill(int color, dmatrix_t lightSource){
//on my computer at home
}
void Polygon::drawLines(int color){
for (int i=0; i<vertexList.size()-1; i++){
	Bresenham(vertexList[i].[1][1], vertexList[i].[2][1], vertexList[i].[1][1], vertexList[i].[2][1], disp, window, GC);
}
//connect last point with first point
Bresenham(vertexList[vertexList.size()-1].[1][1], vertexList[vertexList.size()-1].[2][1], vertexList[0].[1][1], vertexList[0].[2][1], disp, window, GC);
}
/*
taken from my A1 for which i got full marks except for forgetting to comment,
given the x,y of two points, aswell as information about the drawing window, draw that line to the screen.
*/
void Bresenham(int x1, int y1, int x2, int y2, Display *disp, Window window, GC gc){
    int dx, dy, sx, sy;
 	float slope, pitch;
 	//if it is just a single pixel
 	if (x1==x2 && y1==y2){
 		XDrawPoint(disp, window, gc, x1, y1);
 		XFlush(disp);
 		return;
	 }
	 //delta x
	 dx = x2 - x1;
	 //x is decreasing over the line seg
	 if (dx<0){
	 	sx=-1;
	 //x is increasing over the line seg
	 }else{
	 	sx=1;
	 }
	 //delta y
	 dy = y2 - y1;
	 //y is decreasing over the line seg
	 if (dy<0){
	 	sy=-1;
	 }else{//y is increasing over the line seg
	 	sy=1;
	 }
	 //is slope dy/dx or dx/dy
	 if (abs(dy) < abs(dx)){
	 	slope = (float)dy/(float)dx;
	 	pitch = y1 - slope*x1;
	 	//draw the points along the line
	 while(x1!=x2){
	 	XDrawPoint(disp, window, gc, x1, (int)round(slope*x1+pitch));
	 	XFlush(disp);
	 	x1+=sx;
	 }
	 }else{
	 	//same thing but inverted
	 	slope = (float)dx/(float)dy;
	 	pitch = x1 - slope*y1;
	 while (y1!=y2){
	 	XDrawPoint(disp, window, gc, (int)round(slope*y1+pitch), y1);
	 	XFlush(disp);
	 	y1+=sy;
	 }
 }
}
