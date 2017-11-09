/*
Author: Will Jackson
Purpose: This is the header file for the class file that stores a polygon
Date: 11/8/2017
*/
#ifndef POLYGON_H
#define POLYGON_H
#include <list>
#include "matrix.h"
#include <math.h>

class Polygon
{
private:
	list<dmatrix_t> vertexList;
	list<dmatrix_t> normalList;
	dmatrix_t surfaceNormal;

public:
	Polygon(list<dmatrix_t> vertexList);
	~Polygon();
	dmatrix_t getNormal();
	void fill(int color, dmatrix_t lightSource);
	void Bresenham(int x1, int y1, int x2, int y2, Display *disp, Window window, GC gc);
	void drawLines(int color);
};
#endif