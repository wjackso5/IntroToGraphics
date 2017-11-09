/*
Author:Will Jackson
File: this class stores meshes
Date:10/4/17
*/
#include "Mesh.h"


//constructor
Mesh::Mesh(vector<dmatrix_t> vertexList, vector<Polygon> faceList)
{
	vector<dmatrix_t> vertexList = vertexList;
	vector<Polygon> faceList = faceList;
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
Mesh::getVertexList(){
	return vertexList;
}
Mesh::getFaceList(){
	return faceList;
}
