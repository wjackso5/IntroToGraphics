/*
Author:Will Jackson
File: header for meshes
Date:10/4/17
*/
#ifndef MESH_H
#define MESH_H
#include <vector>
#include "matrix.h"
#include <math.h>
#include "Polygon.h"

class Mesh
{
private:
	vector<dmatrix_t> vertexList;
	vector<Polygon> faceList;

public:
	Mesh(vector<dmatrix_t> vertexList, vector<Polygon> faceList );
	~Mesh();
	vector<polygon> getFaceList();
	vector<dmatrix_t> getVertexList();
};
#endif