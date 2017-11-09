/* CS3388 A2 By: William Jackson        
 PURPOSE : Building the synthetic camera for 3D viewing. Display a sphere wire-mesh and a taurus wire-mesh. 
 PREREQUISITES : matrix.h
 This is a modified version of the camera.c file given to us to complete the assignmet, i added code to compute the points of the object and to draw them using x11
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <X11/Xlib.h>
#include "matrix.h"

/*define the radius' for the objects to be created*/
#define SPHERE_R 5.0
#define TORUS_MR 10.0
#define TORUS_R 5.0 



#ifndef NIL
    #define NIL 0
#endif

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/*define the deltas for theta and row used to create the objects*/
#define DELTA_THETA M_PI/360;
#define DELTA_ROW M_PI/180;

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
taken from my A1 for which i got full marks except for forgetting to comment,
given the x,y of two points, aswell as information about the drawing window, draw that line to the screen.
*/
void Bresenham(int x1, int y1, int x2, int y2, Display *disp, Window window, GC gc){
 	printf("x1:%i\n", x1);
    printf("x2:%i\n", x2);
    printf("y1:%i\n", y1);
    printf("y2:%i\n", y2);
    int dx, dy, sx, sy;
 	float slope, pitch;
 	//if it is just a single pixel
 	if (x1==x2 && y1==y2){
 		XDrawPoint(disp, window, gc, x1, y1);
 		XFlush(disp);
 		return;
	 }
	 //delta x
	 dx = x2 - x1;
	 //x is decreasing over the line seg
	 if (dx<0){
	 	sx=-1;
	 //x is increasing over the line seg
	 }else{
	 	sx=1;
	 }
	 //delta y
	 dy = y2 - y1;
	 //y is decreasing over the line seg
	 if (dy<0){
	 	sy=-1;
	 }else{//y is increasing over the line seg
	 	sy=1;
	 }
	 //is slope dy/dx or dx/dy
	 if (abs(dy) < abs(dx)){
	 	slope = (float)dy/(float)dx;
	 	pitch = y1 - slope*x1;
	 	//draw the points along the line
	 while(x1!=x2){
	 	XDrawPoint(disp, window, gc, x1, (int)round(slope*x1+pitch));
	 	XFlush(disp);
	 	x1+=sx;
	 }
	 }else{
	 	//same thing but inverted
	 	slope = (float)dx/(float)dy;
	 	pitch = x1 - slope*y1;
	 while (y1!=y2){
	 	XDrawPoint(disp, window, gc, (int)round(slope*y1+pitch), y1);
	 	XFlush(disp);
	 	y1+=sy;
	 }
 }
}
int main() {
    /*code taken from my assignment 1 now starts,
    I use this to create a window in x11 and draw lines using my implementation of bresenham's algorithm.*/ 
    Display *disp = XOpenDisplay(NIL);
    if (disp) {
        //get black and white colors
        int white = WhitePixel(disp, DefaultScreen(disp));
        int black = BlackPixel(disp, DefaultScreen(disp));
        //Create a simple 512x512 window
        Window window = XCreateSimpleWindow(disp, DefaultRootWindow(disp),0,0,512,512,0, white, white);
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
        /*set the current value of theta and row*/
        double t=0;
        double r=0;
        
        dmatrix_t P1 ; /*set the starting point to draw the sphere. at the very top.*/
       
        dmat_alloc(&P1,4,1) ;

        P1.m[1][1] = cos(t)*sin(r)*SPHERE_R ;
        P1.m[2][1] = sin(t)*sin(r)*SPHERE_R ;
        P1.m[3][1] = cos(r)*SPHERE_R ;
        P1.m[4][1] = 1.0 ;

        dmatrix_t P2 ; /*set the next point to draw the sphere.*/
       
        dmat_alloc(&P2,4,1) ;

        P2.m[1][1] = cos(t)*sin(r)*SPHERE_R ;
        P2.m[2][1] = sin(t)*sin(r)*SPHERE_R ;
        P2.m[3][1] = cos(r)*SPHERE_R ;
        P2.m[4][1] = 1.0 ;
        /*counter var to track whether we should increment theta or row*/
        int c=0;
        /*loop through all the possible values for theta and row to draw our sphere*/
        while (t!=2*M_PI && r!=M_PI){
            /*increment either theta or row*/
            if (c==0&&t!=2*M_PI){
                t = t + DELTA_THETA;
                c = 1; //switch to other case
            }
            else if (c==1&&r!=M_PI){
                r = r + DELTA_ROW;
                c = 0; //switch to other case
            }
            /*set up P2 using new theta/row*/
            P2.m[1][1] = cos(t)*sin(r)*SPHERE_R ;
            P2.m[2][1] = sin(t)*sin(r)*SPHERE_R ;
            P2.m[3][1] = cos(r)*SPHERE_R ;
            P2.m[4][1] = 1.0 ;
            /*multiply the points by the camera matrix and then apply the perspective projections onto those points*/
            P1 = *perspective_projection(dmat_mult(&C,&P1)); 
            P2 = *perspective_projection(dmat_mult(&C,&P2)); 
            /*use bressenhams line drawing algorithm to draw a line from p1 to p2, forming our sphere*/
            Bresenham(round(P1.m[1][1]),round(P1.m[2][1]),round(P2.m[1][1]),round(P2.m[2][1]),disp,window,gc); //This function draws points
            /*set P1 to P2 since we are going to do it all again for the next point*/
            P1 = P2;
        }
        //handle exit input
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
