/*
Author: Will Jackson
Purpose: This is the header file for the class file that stores a polygon
Date: 11/8/2017
*/
#ifndef POLYGON_H
#define POLYGON_H
class Polygon
{
private:
	std::vector<dmatrix_t> vertexList;
	dmatrix_t surfaceNormal;

public:
	Polygon(std::vector<dmatrix_t> vertexList1);
	~Polygon();
	dmatrix_t getNormal();
	std::vector<dmatrix_t> getVertexList();
};
#endif