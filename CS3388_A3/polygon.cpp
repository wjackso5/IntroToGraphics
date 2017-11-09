/*
Author: Will Jackson
Purpose: This is the implementation file for the class file that stores a polygon
Date: 11/8/2017
*/
#include "Polygon.h"
//constructor for a polygon
//gets the list of vertices and then calculates the normals
Polygon::Polygon(std::vector<dmatrix_t> vertexList){
	vertexList=vertexList;
	surfaceNormal = *dcross_product(*mat_sub((vertexList[1]), (vertexList[0]) , *mat_sub(vertexList[2], vertexList[1]));
}
//deconstuctor
Polygon::~Polygon(){
}
//getter for surface normal
dmatrix_t Polygon::getNormal(){
	return surfaceNormal;
}
std::vector<dmatrix_t> Polygon::getVertexList(){
	return vertexList;
}


