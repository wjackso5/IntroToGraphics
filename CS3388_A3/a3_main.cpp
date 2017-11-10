/* CS3388 A2 By: William Jackson   on: oct. 18th, 2017     
 PURPOSE : Building the synthetic camera for 3D viewing. Display a sphere wire-mesh and a torus wire-mesh. 
 PREREQUISITES : matrix.h
 This is a modified version of the camera.c file given to us to complete the assignmet, i added code to compute the points of the object and to draw them using x11
 */
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <X11/Xlib.h>
#include <vector>
#include "Drawer.h"
#include "Mesh.h"
/*define the radius' for the objects to be created
TORUS_RA is the distance from the centre of the torus to the center of the tube
and TORUS_RB is the radius of the tube*/
#define SPHERE_R 4.0
#define TORUS_RA 15.0
#define TORUS_RB 2.5 



#ifndef NIL
    #define NIL 0
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/*define the deltas for theta and row used to create the objects*/
#define DELTA_THETA 15*M_PI/360
#define DELTA_ROW M_PI/180

#define Ex 15.0
#define Ey 15.0
#define Ez 15.0 

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0
#define ASPECT 1.0

#define W  512
#define H  512

void Bresenham(int, int, int, int, Display*, Window, GC);
void CreateSphere(dmatrix_t C);
void CreateTorus(dmatrix_t C);

dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {
    
    dmatrix_t N ; /* Viewing axis */
    
    N = *dmat_normalize(dmat_sub(E,G)) ;
    N.l = 3 ;

    dmatrix_t UP ;
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;
    
    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;
    
    dmatrix_t U ;
    
    U = *dmat_normalize(dcross_product(&UP,&N)) ;
    
    dmatrix_t V ;
    V = *dcross_product(&N,&U) ;
    
    dmatrix_t Mv ; /* Build matrix M_v */
    dmat_alloc(&Mv,4,4) ;
    
    Mv.m[1][1] = U.m[1][1] ; 
    Mv.m[1][2] = U.m[2][1] ; 
    Mv.m[1][3] = U.m[3][1] ; 
    Mv.m[1][4] = -1.0*((*E).m[1][1]*U.m[1][1] + (*E).m[2][1]*U.m[2][1] + (*E).m[3][1]*U.m[3][1]) ;
    
    Mv.m[2][1] = V.m[1][1] ; 
    Mv.m[2][2] = V.m[2][1] ;    
    Mv.m[2][3] = V.m[3][1] ; 
    Mv.m[2][4] = -1.0*((*E).m[1][1]*V.m[1][1] + (*E).m[2][1]*V.m[2][1] + (*E).m[3][1]*V.m[3][1]) ;
    
    Mv.m[3][1] = N.m[1][1] ; 
    Mv.m[3][2] = N.m[2][1] ; 
    Mv.m[3][3] = N.m[3][1] ; 
    Mv.m[3][4] = -1.0*((*E).m[1][1]*N.m[1][1] + (*E).m[2][1]*N.m[2][1] + (*E).m[3][1]*N.m[3][1]) ;
    
    Mv.m[4][1] = 0.0 ; 
    Mv.m[4][2] = 0.0 ; 
    Mv.m[4][3] = 0.0 ; 
    Mv.m[4][4] = 1.0 ;
    
    dmatrix_t Mp ; /* Build matrix Mp */
    dmat_alloc(&Mp,4,4) ;
    Mp = *dmat_identity(&Mp) ;
    
    float a = -1.0*(FP + NP)/(FP - NP) ;
    float b = -2.0*(FP*NP)/(FP - NP) ;
    
    Mp.m[1][1] = NP ;
    Mp.m[2][2] = NP ;
    Mp.m[3][3] = a ;
    Mp.m[3][4] = b ;
    Mp.m[4][3] = -1.0 ;
    Mp.m[4][4] = 0.0 ;
    
    /* Build matrices T_1 and S_1 */
    
    /* Work out coordinates of near plane corners */
    
    float top = NP*tan(M_PI/180.0*THETA/2.0) ;
    float right = ASPECT*top ;
    float bottom = -top ;
    float left = -right ;
   
    dmatrix_t T1 ;
    dmat_alloc(&T1,4,4) ;
    
    T1 = *dmat_identity(&T1) ;
    T1.m[1][4] = -(right + left)/2.0 ;
    T1.m[2][4] = -(top + bottom)/2.0 ;

    dmatrix_t S1 ;
    dmat_alloc(&S1,4,4) ;
    
    S1 = *dmat_identity(&S1) ;
    S1.m[1][1] = 2.0/(right - left) ;
    S1.m[2][2] = 2.0/(top - bottom) ;

    /* Build matrices T2, S2, and W2 */
    
    dmatrix_t T2 ;
    dmatrix_t S2 ;
    dmatrix_t W2 ;
    
    dmat_alloc(&T2,4,4) ;
    dmat_alloc(&S2,4,4) ;
    dmat_alloc(&W2,4,4) ;
    
    T2 = *dmat_identity(&T2) ;
    S2 = *dmat_identity(&S2) ;
    W2 = *dmat_identity(&W2) ;
    
    T2.m[1][4] = 1.0 ;
    T2.m[2][4] = 1.0 ;

    S2.m[1][1] = W/2.0 ;
    S2.m[2][2] = H/2.0 ;
    
    W2.m[2][2] = -1.0 ;
    W2.m[2][4] = (double)H ;
    
    return dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,dmat_mult(&S1,dmat_mult(&T1,dmat_mult(&Mp,&Mv)))))) ;
}

dmatrix_t *perspective_projection(dmatrix_t *P) {

    (*P).m[1][1] /= (*P).m[4][1] ;
    (*P).m[2][1] /= (*P).m[4][1] ;
    (*P).m[3][1] /= (*P).m[4][1] ;
    (*P).m[4][1] /= (*P).m[4][1] ;

    return P ;
}
/*
This function draws a sphere using our implementation of bresenham by iterating through all theta between 0...2*Pi
and through all row between 0...Pi. It takes four parameters, the first three are all information about the X11 display so
that we can call our implementation of bresenham. The fourth is the camera Matrix which we need to preform perspective projections.
The radius of the sphere is declared in the header of this file as a constant as well as the deltas for row and theta.
*/
Mesh createSphere(dmatrix_t C){
    //declare Vector and face List for the mesh to be made
    std::vector<dmatrix_t> vertexList;
    std::vector<Polygon> faceList;
    /*set the current value of theta and row*/
        double t=0;
        double r=0;
        
        dmatrix_t P0 ; /*set the starting point to draw the sphere. at the very top.*/
        dmat_alloc(&P0,4,1) ;
        dmatrix_t P1 ; /*set the starting point to draw the sphere. at the very top.*/
        dmat_alloc(&P1,4,1) ;
        dmatrix_t P2 ; /*set the starting point to draw the sphere. at the very top.*/
        dmat_alloc(&P2,4,1) ;
        dmatrix_t P3 ; /*set the starting point to draw the sphere. at the very top.*/
        dmat_alloc(&P3,4,1) ;

        /*loop through all the possible values for theta and row to draw our sphere*/
        while (t < 2*M_PI){
            t = t + DELTA_THETA;
            r=0;
            while (r < M_PI){
            r = r + DELTA_ROW;
            /*set up P0 using new theta/row*/
            P0.m[1][1] = cos(t)*sin(r)*SPHERE_R;
            P0.m[2][1] = sin(t)*sin(r)*SPHERE_R;
            P0.m[3][1] = cos(r)*SPHERE_R;
            P0.m[4][1] = 1.0 ;
            // for P1
            P1.m[1][1] = cos(t+DELTA_THETA)*sin(r)*SPHERE_R;
            P1.m[2][1] = sin(t+DELTA_THETA)*sin(r)*SPHERE_R;
            P1.m[3][1] = cos(r)*SPHERE_R;
            P1.m[4][1] = 1.0 ;
            // for P2
            P2.m[1][1] = cos(t+DELTA_THETA)*sin(r+DELTA_ROW)*SPHERE_R;
            P2.m[2][1] = sin(t+DELTA_THETA)*sin(r+DELTA_ROW)*SPHERE_R;
            P2.m[3][1] = cos(r+DELTA_ROW)*SPHERE_R;
            P2.m[4][1] = 1.0 ;
            // for P3
            P3.m[1][1] = cos(t)*sin(r+DELTA_ROW)*SPHERE_R;
            P3.m[2][1] = sin(t)*sin(r+DELTA_ROW)*SPHERE_R;
            P3.m[3][1] = cos(r+DELTA_ROW)*SPHERE_R;
            P3.m[4][1] = 1.0 ;
            //apply perspective projections
            dmatrix_t Pp0; 
            dmat_alloc(&Pp0,4,1);
            dmatrix_t Pp1; 
            dmat_alloc(&Pp1,4,1);
            dmatrix_t Pp2; 
            dmat_alloc(&Pp2,4,1);
            dmatrix_t Pp3; 
            dmat_alloc(&Pp3,4,1);
            Pp0 = *perspective_projection(dmat_mult(&C,&P0)); 
            Pp1 = *perspective_projection(dmat_mult(&C,&P1));
            Pp2 = *perspective_projection(dmat_mult(&C,&P2)); 
            Pp3 = *perspective_projection(dmat_mult(&C,&P3));
            //add Pp0 to vertex list.
            vertexList.push_back(Pp0);
            //poly list for this polygon
            std::vector<dmatrix_t> polyVertexList;
            polyVertexList.push_back(Pp0);
            polyVertexList.push_back(Pp1);
            polyVertexList.push_back(Pp2);
            polyVertexList.push_back(Pp3);
            //create the polygon
            Polygon poly = Polygon(polyVertexList);
            //add polygon to faceList
            faceList.push_back(poly);
            }
        }
        //now that we have populated the lists, we will make the mesh for the sphere
        Mesh sphere = Mesh(vertexList, faceList);
        return sphere;
}
/*
This function creates a torus using our by iterating through all theta between 0...2*Pi
*/
Mesh createTorus(dmatrix_t C){
    /*set the current value of theta and row*/
        std::vector<dmatrix_t> vertexList;
        std::vector<Polygon> faceList;
        double t=0;
        double r=0;
        //declare points for polygon
        dmatrix_t P0 ;
        dmat_alloc(&P0,4,1) ;
        dmatrix_t P1 ;
        dmat_alloc(&P1,4,1) ;
        dmatrix_t P2 ;
        dmat_alloc(&P2,4,1) ;
        dmatrix_t P3 ; 
        dmat_alloc(&P3,4,1) ;

        /*loop through all the possible values for theta and row to genertae our torus*/
        while (t < 2*M_PI){
            t = t + DELTA_THETA;
            r=0;
            while (r < 2*M_PI){
            
            r = r + DELTA_ROW;
            //points on the polygon
            P0.m[1][1] = (TORUS_RA+(TORUS_RB*cos(t)))*cos(r);
            P0.m[2][1] = (TORUS_RA+(TORUS_RB*cos(t)))*sin(r);
            P0.m[3][1] = TORUS_RB*sin(t);
            P0.m[4][1] = 1.0 ;
            // for P1
            P1.m[1][1] = (TORUS_RA+(TORUS_RB*cos(t+DELTA_THETA)))*cos(r);
            P1.m[2][1] = (TORUS_RA+(TORUS_RB*cos(t+DELTA_THETA)))*sin(r);
            P1.m[3][1] = TORUS_RB*sin(t+DELTA_THETA);
            P1.m[4][1] = 1.0 ;
            // for P2
            P2.m[1][1] = (TORUS_RA+(TORUS_RB*cos(t+DELTA_THETA)))*cos(r+DELTA_ROW);
            P2.m[2][1] = (TORUS_RA+(TORUS_RB*cos(t+DELTA_THETA)))*sin(r+DELTA_ROW);
            P2.m[3][1] = TORUS_RB*sin(t+DELTA_THETA);
            P2.m[4][1] = 1.0 ;
            // for P3
            P3.m[1][1] = (TORUS_RA+(TORUS_RB*cos(t)))*cos(r+DELTA_ROW);
            P3.m[2][1] = (TORUS_RA+(TORUS_RB*cos(t)))*sin(r+DELTA_ROW);
            P3.m[3][1] = TORUS_RB*sin(t);
            P3.m[4][1] = 1.0 ;
            //apply perspective projections
            dmatrix_t Pp0; 
            dmat_alloc(&Pp0,4,1);
            dmatrix_t Pp1; 
            dmat_alloc(&Pp1,4,1);
            dmatrix_t Pp2; 
            dmat_alloc(&Pp2,4,1);
            dmatrix_t Pp3; 
            dmat_alloc(&Pp3,4,1);
            Pp0 = *perspective_projection(dmat_mult(&C,&P0)); 
            Pp1 = *perspective_projection(dmat_mult(&C,&P1));
            Pp2 = *perspective_projection(dmat_mult(&C,&P2)); 
            Pp3 = *perspective_projection(dmat_mult(&C,&P3));
            //add Pp0 to vertex list.
            vertexList.push_back(Pp0);
            //poly list for this polygon
            std::vector<dmatrix_t> polyVertexList;
            polyVertexList.push_back(Pp0);
            polyVertexList.push_back(Pp1);
            polyVertexList.push_back(Pp2);
            polyVertexList.push_back(Pp3);
            //create the polygon
            Polygon poly = Polygon(polyVertexList);
            //add polygon to faceList
            faceList.push_back(poly);
        }
        //now that we have populated the lists, we will make the mesh for the torus
        Mesh torus = Mesh(vertexList, faceList);
        return torus;
    }



}
/*
use code from A1 to set up an X11 window then sets up the camera matrix and then draws a sphere and a torus onto the window, waits for input q
and on that input closes. 
*/
int main() {
    /*code taken from my assignment 1 now starts,
    I use this to create a window in x11 and draw lines using my implementation of bresenham's algorithm.*/ 
    Display *disp = XOpenDisplay(NIL);
    if (disp) {
        //get black and white colors
        int white = WhitePixel(disp, DefaultScreen(disp));
        int black = BlackPixel(disp, DefaultScreen(disp));
        //Create a simple 512x512 window
        Window window = XCreateSimpleWindow(disp, DefaultRootWindow(disp),0,0,W,H,0, white, white);
        XSelectInput(disp,window,StructureNotifyMask);
        //Make the window appear
        XMapWindow(disp, window);
        //graphics context
        GC gc = XCreateGC(disp, window, 0, NIL);
        //set draw to black
        XSetForeground(disp, gc, black);
        //Wait to get events.
        for(;;){
            XEvent ev;
            XNextEvent(disp, &ev);
            if (ev.type == MapNotify){
            break;
            }
        }
        /*display window is now ready*/

        dmatrix_t E ; /* The centre of projection for the camera */
        
        dmat_alloc(&E,4,1) ;
        
        E.m[1][1] = Ex ;
        E.m[2][1] = Ey ;
        E.m[3][1] = Ez ;
        E.m[4][1] = 1.0 ;
        
        dmatrix_t G ; /* Point gazed at by camera */
        
        dmat_alloc(&G,4,1) ;
        
        G.m[1][1] = Gx ;
        G.m[2][1] = Gy ;
        G.m[3][1] = Gz ;
        G.m[4][1] = 1.0 ;

        dmatrix_t C ; /* The camera matrix */

        dmat_alloc(&C,4,4);
        C = *build_camera_matrix(&E,&G) ;
        Mesh sphere = createSphere(C);
        Drawer drawer = Drawer(disp, window, gc, C, C);
        drawer.drawMesh(sphere);
        
        char inp='0';
        while(inp!='q'){
           printf("press q+enter to exit: ");
           inp = getchar();
           if (inp!='q'){
           printf("invalid input\n");
           }
        }	
	
    }
    
    return 0;
}
