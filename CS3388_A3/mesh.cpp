/*
Author:Will Jackson
File: this class stores meshes
Date:10/4/17
*/
#include "Mesh.h"


//constructor
Mesh::Mesh(int drawlines, int polyfill, int hideback, vector<dmatrix> vertexList, vector<polygon> faceList, dmatrix camera, dmatrix lightSource, int color)
{
	vector<dmatrix> vertexList = vertexList;
	vector<polygon> faceList = faceList;
	int drawlines, polyfill, hideback;
	//show the mesh based on option params
	for (int i=0; i<faceList.size(); i++;){
		if ((hideback==1) && (/*angle betweeen face.getNormal() and camera matrix is gt 90*/)){
		//do nothing as the surface is hidden
		}else{
			if (polyfill=1){faceList[i].fill(color,lightSource);}
			if (drawlines=1){faceList[i].drawlines(color);}
		}
	}
}
//deconstructor
Mesh::~Mesh()
{
}
