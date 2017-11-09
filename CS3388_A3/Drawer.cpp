/*Author: Will Jackson
Purpose: This is the header file for the drawer class takes care of the visualization of our meshes.
Date: 11/8/2017
*/
#include "Drawer.h"


//constructor
Drawer::Drawer(Display *disp1, Window window1, GC gc1, dmatrix_t camera1, dmatrix_t light1)
{
	disp=disp1;
	window=window1;
	gc=gc1;
	camera = camera1;
	light = light1;
}
//deconstructor
Drawer::~Drawer()
{
}
void fillPolygon(Polygon p){

}

/*
taken from my A1 for which i got full marks except for forgetting to comment,
given the x,y of two points, aswell as information about the drawing window, draw that line to the screen.
*/
void Bresenham(int x1, int y1, int x2, int y2){

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

//draws the lines of the polygon
void drawPolygon(Polygon p){
for (int i=0; i<p.getVertexList().size()-1; i++){
	Bresenham(p.getVertexList().at(i).m[1][1], p.getVertexList().at(i).m[2][1], p.getVertexList().at(i).m[1][1], p.getVertexList().at(i).m[2][1]);
}
//draw last line between final point and first point
//connect last point with first point
Bresenham(p.getVertexList().at(p.getVertexList().size()-1).m[1][1], p.getVertexList().at(p.getVertexList().size()-1).m[2][1], p.getVertexList().at(0)[1][1], p.getVertexList().at(0).m[2][1]);
}