/*
Author: Will Jackson
Purpose: This is the header file for the class file that stores a polygon
Date: 11/8/2017
*/
#ifndef POLYGON_H
#define POLYGON_H
#include <list>
#include "matrix.h"
class Polygon
{
private:
	list<dmatrix_t> vertexList;
	list<dmatrix_t> normalList;
	dmatrix_t surfaceNormal;
public:
	Polygon(list<dmatrix_t> vertexList);
	void fill(int color, dmatrix_t lightSource);
	void drawLines(int color);
};
#endif