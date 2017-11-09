/*Author: Will Jackson
Purpose: This is the header file for the drawer class takes care of the visualization of our meshes.
Date: 11/8/2017
*/
#ifndef DRAWER_H
#define DRAWER_H
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <X11/Xlib.h>
#include <vector>
#include "Mesh.h"
extern "C" {
#include "matrix.h"
}

class Drawer
{
private:
Display *disp;
Window window;
GC gc;
dmatrix_t camera, light;
public:
	Drawer(Display *disp1, Window window1, GC gc1, dmatrix_t camera1, dmatrix_t light1);
	~Drawer();

	void fillPolygon(Polygon p);
	void Bresenham(int x1, int y1, int x2, int y2);
	void drawPolygon(Polygon p);
};
#endif