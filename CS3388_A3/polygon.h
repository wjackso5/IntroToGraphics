/*
Author: Will Jackson
Purpose: This is the header file for the class file that stores a polygon
Date: 11/8/2017
*/
#ifndef POLYGON_H
#define POLYGON_H
#include <vector>
#include "matrix.h"
#include <math.h>

class Polygon
{
private:
	vector<dmatrix_t> vertexList;
	dmatrix_t surfaceNormal;

public:
	Polygon(vector<dmatrix_t> vertexList);
	~Polygon();
	dmatrix_t getNormal();
	vector<dmatrix_t> getVertexList();
};
#endif