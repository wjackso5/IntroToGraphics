/*
Author:Will Jackson
File: header for meshes
Date:10/4/17
*/
#ifndef MESH_H
#define MESH_H

#include "Polygon.h"


class Mesh
{
private:
	std::vector<dmatrix_t> vertexList;
	std::vector<Polygon> faceList;

public:
	Mesh(std::vector<dmatrix_t> vertexList1, std::vector<Polygon> faceList1);
	~Mesh();
	std::vector<Polygon> getFaceList();
	std::vector<dmatrix_t> getVertexList();
};
#endif