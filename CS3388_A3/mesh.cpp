/*
Author:Will Jackson
File: this class stores meshes
Date:10/4/17
*/
#include "Mesh.h"


//constructor
Mesh::Mesh(std::vector<dmatrix_t> vertexList1, std::vector<Polygon> faceList1)
{
	std::vector<dmatrix_t> vertexList = vertexList1;
	std::vector<Polygon> faceList = faceList1;
	/*this needs to go in drawer or main
	//show the mesh based on option params
	for (int i=0; i<faceList.size(); i++;){
		if ((hideback==1) && (angle betweeen face.getNormal() and camera matrix is gt 90)){
		//do nothing as the surface is hidden
		}else{
			if (polyfill=1){faceList[i].fill(color,lightSource);}
			if (drawlines=1){faceList[i].drawlines(color);}
		}
	}*/

}
//deconstructor
Mesh::~Mesh()
{
}
std::vector<dmatrix_t> Mesh::getVertexList(){
	return vertexList;
}
std::vector<Polygon> Mesh::getFaceList(){

	return faceList;
}
