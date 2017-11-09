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
extern "C" {
#include "matrix.h"
}

class Mesh
{
private:
	std::vector<dmatrix_t> vertexList;
	std::vector<Polygon> faceList;

public:
	Mesh(std::vector<dmatrix_t> vertexList, std::vector<Polygon> faceList );
	~Mesh();
	std::vector<polygon> getFaceList();
	std::vector<dmatrix_t> getVertexList();
};
#endif