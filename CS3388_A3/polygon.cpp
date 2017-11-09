/*
Author: Will Jackson
Purpose: This is the implementation file for the class file that stores a polygon
Date: 11/8/2017
*/
#include "polygon.h"
#include <math.h>
#include "matrix.h"
//constructor for a polygon
//gets the list of vertices and then calculates the normals
Polygon(list<dmatrix_t> vertexList){
	vertexList=vertexList;
	surfaceNormal = *dcross_product(*mat_sub(*vertexList[1], *vertexList[0]) , *mat_sub(*vertexList[2], *vertexList[1]))
}
void fill(int color, dmatrix_t lightSource){

}
void drawLines(int color){

}

#endif