/*
Author:Will Jackson
File: header for meshes
Date:10/4/17
*/
#ifndef MESH_H
#define MESH_H
#include <list>
#include "matrix.h"
#include <math.h>

class Mesh
{
private:
	vector<dmatrix> vertexList;
	vector<polygon> faceList;

public:
	Polygon(list<dmatrix_t> vertexList);
	~Polygon();
	dmatrix_t getNormal();
	void getVertexList();
	void getFaceList();
	void fill(int color, dmatrix_t lightSource);
	void Bresenham(int x1, int y1, int x2, int y2, Display *disp, Window window, GC gc);
	void drawLines(int color);
};
#endif