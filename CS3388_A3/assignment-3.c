/* Program which uses X11 and implements the synthetic camera and traces parametric functions also fills the shapes with shaded polygons.
based off of assignment-2-solution.c*/

#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#define FALSE 0
#define TRUE  1

#define Ex 20.0
#define Ey 20.0
#define Ez 100.0

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 10.0
#define FP 50.0

#define THETA 90.0

#define W  1400
#define H  800

#define POSX  0
#define POSY  0

#define ASPECT (double)W/(double)H

#define TRIANGLE 3

#define EPSILON 0.0001

#define M_PI 3.1415926535897932384626433832795028

Display *InitX(Display *d,Window *w,int *s) {
    
  d = XOpenDisplay(NULL) ;
  if(d == NULL) {
    printf("Cannot open display\n") ;
    exit(1) ;
  }
  *s = DefaultScreen(d) ;
  *w = XCreateSimpleWindow(d,RootWindow(d, *s),POSX,POSY,W,H,1,BlackPixel(d, *s), WhitePixel(d, *s)) ;
  Atom delWindow = XInternAtom(d,"WM_DELETE_WINDOW",0) ;
  XSetWMProtocols(d,*w,&delWindow,1) ;
  XSelectInput(d,*w,ExposureMask | KeyPressMask) ;
  XMapWindow(d,*w) ;
  return(d) ;
}


void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {
    
  XSetForeground(d,*gc,r << 16 | g << 8 | b) ;
}


void SetPixelX(Display *d, Window w, int s, int i, int j) {
    
  XDrawPoint(d,w,DefaultGC(d,s),i,j) ;
}


void QuitX(Display *d, Window w) {
    
  XDestroyWindow(d,w) ;
  XCloseDisplay(d) ;
}


void exchangeInt(int *a, int *b)

{ int t ;
    
  t = *a ;
  *a = *b ;
  *b = t ;
}

void bresenham(Display *d, Window w, int s, int x1, int y1, int x2, int y2)

{ int Transx, Transy, Pi, Dx, Dy, Two_Dx, Two_Dy, i, Inc1stcoord, Inc2ndcoord, Exchange ;
  
  Exchange = FALSE ;
  Inc1stcoord = 1 ;
  Inc2ndcoord = 1 ;
    
  Transx = -x1 ;
  Transy = -y1 ;
    
  x1 = 0 ;
  y1 = 0 ;
    
  x2 += Transx ;
  y2 += Transy ;
    
  Dx = x2 ;
  Dy = y2 ;
    
  Two_Dx = 2*x2 ;
  Two_Dy = 2*y2 ;
    
  if (Dy < 0) {
    Inc2ndcoord = -1 ;
    Dy *= -1 ;
    Two_Dy *= -1 ;
  }
  if (Dx < 0) {
    Inc1stcoord = -1 ;
    Dx *= -1 ;
    Two_Dx *= -1 ;
  }
  if (Dy > Dx) {
    Exchange = TRUE ;
    exchangeInt(&Two_Dx,&Two_Dy) ;
    exchangeInt(&Dx,&Dy) ;
    exchangeInt(&Inc1stcoord,&Inc2ndcoord) ;
  }
  Pi = Two_Dy - Dx ;
  if (Exchange) {
    SetPixelX(d,w,s,y1 - Transx,x1 - Transy) ;
  }
  else {
    SetPixelX(d,w,s,x1 - Transx,y1 - Transy) ;
  }
  for (i = 0 ; i < Dx ; i++) {
    if (Pi < 0) {
      Pi += Two_Dy ;
    }
    else {
      Pi += Two_Dy - Two_Dx ;
      y1 += Inc2ndcoord ;
    }
    x1 += Inc1stcoord ;
    if (Exchange) {
      SetPixelX(d,w,s,y1 - Transx,x1 - Transy) ;
    }
    else {
      SetPixelX(d,w,s,x1 - Transx,y1 - Transy) ;
    }
  }
}
//implementation of line scan polygon filling algorithm
void polyfill(Display *d, Window w, int s, int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4){
	//find ymin
	double i1, i2, m1, m2, m3, m4;
	int i1f=0;
	int i2f=0;
	int ymin=y1;
	if (y2<ymin)ymin=y2;
	if (y3<ymin)ymin=y3;
	if (y4<ymin)ymin=y4;
	//find ymax
	int ymax=y1;
	if (y2<ymax)ymax=y2;
	if (y3<ymax)ymax=y3;
	if (y4<ymax)ymax=y4;
	int scanline = ymax;
	//calculate slopes
	m1=(y2-y1)/(x2-x1);
	m2=(y3-y2)/(x3-x2);
	m3=(y4-y3)/(x4-x3);
	m4=(y1-y4)/(x1-x4);
	while (scanline >= ymin){
		//calculate instersections with scanline and fill
			if  (m1>2147483634||m1<-2147483634){i1=x1;}else{
			i1 = ((scanline-y1)/m1)+x1;}
			if (m3>2147483634||m3<-2147483634){i2=x3;}else{
			i2 = ((scanline-y3)/m3)+x3;}
			
			bresenham(d,w,s,(int)i1,scanline,(int)i2,scanline);
		
		  scanline = scanline -1;
	}
}

dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {
    
  dmatrix_t N, UP, U, V, Mv, Mp, T1, S1, T2, S2, W2 ; 
    
  N = *dmat_normalize(dmat_sub(E,G)) ;
  N.l = 3 ;

  dmat_alloc(&UP,4,1) ;
  UP.l = 3 ;
  UP.m[1][1] = UPx ;
  UP.m[2][1] = UPy ;
  UP.m[3][1] = UPz ;
  UP.m[4][1] = 0.0 ;
    
  U = *dmat_normalize(dcross_product(&UP,&N)) ;
  V = *dcross_product(&N,&U) ;

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
   
  dmat_alloc(&T1,4,4) ;
  T1 = *dmat_identity(&T1) ;
  T1.m[1][4] = -(right + left)/2.0 ;
  T1.m[2][4] = -(top + bottom)/2.0 ;

  dmat_alloc(&S1,4,4) ;
  S1 = *dmat_identity(&S1) ;
  S1.m[1][1] = 2.0/(right - left) ;
  S1.m[2][2] = 2.0/(top - bottom) ;

  /* Build matrices T2, S2, and W2 */
    
  dmat_alloc(&T2,4,4) ;
  dmat_alloc(&S2,4,4) ;
  dmat_alloc(&W2,4,4) ;
    
  T2 = *dmat_identity(&T2) ;
  S2 = *dmat_identity(&S2) ;
  W2 = *dmat_identity(&W2) ;
    
  T2.m[1][4] = 1.0 ;
  T2.m[2][4] = 1.0 ;

  S2.m[1][1] = (double)W/2.0 ;
  S2.m[2][2] = (double)H/2.0 ;
    
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


dmatrix_t *sphere_point(dmatrix_t *P, double r, double theta, double rho) {

  (*P).m[1][1] = r*cos(theta)*sin(rho) ;
  (*P).m[2][1] = r*sin(theta)*sin(rho) ;
  (*P).m[3][1] = r*cos(rho) ;
  (*P).m[4][1] = 1.0 ;

  return P ;
}

dmatrix_t *cone_point(dmatrix_t *P, double r, double theta, double rho) {

  (*P).m[1][1] = r*cos(theta) ;
  (*P).m[2][1] = r*sin(theta) ;
  (*P).m[3][1] = r;
  (*P).m[4][1] = 1.0 ;

  return P ;
}


dmatrix_t *torus_point(dmatrix_t *P, double a, double c, double theta, double rho) {

  (*P).m[1][1] = (c + a*cos(rho))*cos(theta) ;
  (*P).m[2][1] = (c + a*cos(rho))*sin(theta) ;
  (*P).m[3][1] = a*sin(rho) ;
  (*P).m[4][1] = 1.0 ;

  return P ;
}


void shade_wiremesh_sphere(Display *d, Window w, int s, dmatrix_t *C, dmatrix_t *L, double radius, double theta_lower_bound, double theta_upper_bound, double rho_lower_bound, double rho_upper_bound, double delta_theta, double delta_rho) {
  int i ;
  double theta, rho ;
  //points of the current polygon
  dmatrix_t P[TRIANGLE+1] ;
  dmatrix_t N;
  dmat_alloc(&N,4,1)	 ;

  for (i = 0 ; i < TRIANGLE+1 ; i++) {
    dmat_alloc(&P[i],4,1) ;
  }
  
  for (theta = theta_lower_bound ; theta + delta_theta < theta_upper_bound + EPSILON ; theta += delta_theta) {
    for (rho = rho_lower_bound ; rho + delta_rho < rho_upper_bound + EPSILON ; rho += delta_rho) {
    	
    	//points of the polygon
		P[0] = *perspective_projection(dmat_mult(C, sphere_point(&P[0], radius, theta, rho)));
		P[1] = *perspective_projection(dmat_mult(C, sphere_point(&P[1], radius, theta + delta_theta, rho)));
		P[2] = *perspective_projection(dmat_mult(C, sphere_point(&P[2], radius, theta + delta_theta, rho + delta_rho)));
		P[3] = *perspective_projection(dmat_mult(C, sphere_point(&P[3], radius, theta, rho + delta_rho)));
		//calculate surface normal

		N = *dcross_product(from_homogeneous(dmat_sub(&P[1], &P[0])), from_homogeneous(dmat_sub(&P[2], &P[1])));
		if(N.m[1][1]==0&&N.m[2][1]==0&&N.m[3][1]==0){
		N = *dcross_product(from_homogeneous(dmat_sub(&P[1], &P[3])) , from_homogeneous(dmat_sub(&P[2], &P[1])));
		}
		
		//draw outlines and shade using hideback
		if (acos((Ex*N.m[1][1]+Ey*N.m[2][1]+Ez*N.m[3][1])/(sqrt(N.m[1][1]*N.m[1][1]+N.m[2][1]*N.m[2][1]+N.m[3][1]*N.m[3][1])*sqrt((Ez*Ez)+(Ey*Ey)+(Ex*Ex))))<=(M_PI/2)){
		  bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[1].m[1][1],(int)P[1].m[2][1]) ;
      bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[2].m[1][1],(int)P[2].m[2][1]) ;
        //polyfill(d,w,s,(int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[1].m[1][2],P[1].m[2][1],(int)P[2].m[1][1], (int)P[2].m[2][1],(int)P[3].m[1][1], (int)P[3].m[2][1]);
		}	
    }
  }
  for (i = 0 ; i < TRIANGLE+1 ; i++) {
    free_dmatrix(P[i].m,1,P[i].l,1,P[i].c) ;
  }
}

void shade_wiremesh_cone(Display *d, Window w, int s, dmatrix_t *C, dmatrix_t *L, double radius, double theta_lower_bound, double theta_upper_bound, double rho_lower_bound, double rho_upper_bound, double delta_theta, double delta_rho) {
  int i ;
  double theta, rho ;
  //points of the current polygon
  dmatrix_t P[TRIANGLE+1] ;
  dmatrix_t N;
  dmat_alloc(&N,4,1)   ;

  for (i = 0 ; i < TRIANGLE+1 ; i++) {
    dmat_alloc(&P[i],4,1) ;
  }
  
  for (theta = theta_lower_bound ; theta + delta_theta < theta_upper_bound + EPSILON ; theta += delta_theta) {
    for (rho = rho_lower_bound ; rho + delta_rho < rho_upper_bound + EPSILON ; rho += delta_rho) {
      
      //points of the polygon
    P[0] = *perspective_projection(dmat_mult(C, cone_point(&P[0], radius, theta, rho)));
    P[1] = *perspective_projection(dmat_mult(C, cone_point(&P[1], radius, theta + delta_theta, rho)));
    P[2] = *perspective_projection(dmat_mult(C, cone_point(&P[2], radius, theta + delta_theta, rho + delta_rho)));
    P[3] = *perspective_projection(dmat_mult(C, cone_point(&P[3], radius, theta, rho + delta_rho)));
    //calculate surface normal

    N = *dcross_product(from_homogeneous(dmat_sub(&P[1], &P[0])), from_homogeneous(dmat_sub(&P[2], &P[1])));
    if(N.m[1][1]==0&&N.m[2][1]==0&&N.m[3][1]==0){
    N = *dcross_product(from_homogeneous(dmat_sub(&P[1], &P[3])) , from_homogeneous(dmat_sub(&P[2], &P[1])));
    }
    
    //draw outlines and shade using hideback
    if (acos((Ex*N.m[1][1]+Ey*N.m[2][1]+Ez*N.m[3][1])/(sqrt(N.m[1][1]*N.m[1][1]+N.m[2][1]*N.m[2][1]+N.m[3][1]*N.m[3][1])*sqrt((Ez*Ez)+(Ey*Ey)+(Ex*Ex))))<=(M_PI/2)){
      bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[1].m[1][1],(int)P[1].m[2][1]) ;
      bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[2].m[1][1],(int)P[2].m[2][1]) ;
        //polyfill(d,w,s,(int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[1].m[1][2],P[1].m[2][1],(int)P[2].m[1][1], (int)P[2].m[2][1],(int)P[3].m[1][1], (int)P[3].m[2][1]);
    } 
    }
  }
  for (i = 0 ; i < TRIANGLE+1 ; i++) {
    free_dmatrix(P[i].m,1,P[i].l,1,P[i].c) ;
  }
}


void shade_wiremesh_torus(Display *d, Window w, int s, dmatrix_t *C, dmatrix_t *L, double a, double c, double theta_lower_bound, double theta_upper_bound, double rho_lower_bound, double rho_upper_bound, double delta_theta, double delta_rho) {
  int i ;
  double theta, rho ;
  //points of the current polygon
  dmatrix_t P[TRIANGLE+1] ;
  dmatrix_t N;
  dmat_alloc(&N,4,1);

  for (i = 0 ; i < TRIANGLE+1 ; i++) {
    dmat_alloc(&P[i],4,1) ;
  }
  
  for (theta = theta_lower_bound ; theta + delta_theta < theta_upper_bound + EPSILON ; theta += delta_theta) {
    for (rho = rho_lower_bound ; rho + delta_rho < rho_upper_bound + EPSILON ; rho += delta_rho) {
    	
    	//points of the polygon
		P[0] = *perspective_projection(dmat_mult(C,torus_point(&P[0],a,c,theta,rho))) ; 
    P[1] = *perspective_projection(dmat_mult(C,torus_point(&P[1],a,c,theta+delta_theta,rho))) ;
    P[2] = *perspective_projection(dmat_mult(C,torus_point(&P[2],a,c,theta+delta_rho,rho+delta_rho))) ;
		P[3] = *perspective_projection(dmat_mult(C,torus_point(&P[3], a,c, theta, rho + delta_rho)));
		//calculate surface normal

		N = *dcross_product(from_homogeneous(dmat_sub(&P[1], &P[0])), from_homogeneous(dmat_sub(&P[2], &P[1])));
		if(N.m[1][1]==0&&N.m[2][1]==0&&N.m[3][1]==0){
		N = *dcross_product(from_homogeneous(dmat_sub(&P[1], &P[3])) , from_homogeneous(dmat_sub(&P[2], &P[1])));
		}
		
		//draw outlines and shade using hideback
		if (acos((Ex*N.m[1][1]+Ey*N.m[2][1]+Ez*N.m[3][1])/(sqrt(N.m[1][1]*N.m[1][1]+N.m[2][1]*N.m[2][1]+N.m[3][1]*N.m[3][1])*sqrt((Ez*Ez)+(Ey*Ey)+(Ex*Ex))))<=(M_PI/2)){
			//get the correct color shading
			double Lx=Ex+150.0;
			double Ly=Ey+100.0;
			double Lz=Ez;

			double angle = acos((Lx*N.m[1][1]+Lx*N.m[2][1]+Lz*N.m[3][1])/(sqrt(N.m[1][1]*N.m[1][1]+N.m[2][1]*N.m[2][1]+N.m[3][1]*N.m[3][1])*sqrt((Lx*Lx)+(Ly*Ly)+(Lz*Lz))))	;	
			SetCurrentColorX(d,&(DefaultGC(d,s)),0,0,(angle*(2/M_PI)*255));
			bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[1].m[1][1],(int)P[1].m[2][1]) ;
      bresenham(d,w,s,(int)P[0].m[1][1],(int)P[0].m[2][1],(int)P[2].m[1][1],(int)P[2].m[2][1]) ;
        	//polyfill(d,w,s,(int)P[0].m[1][1], (int)P[0].m[2][1], (int)P[1].m[1][2],P[1].m[2][1],(int)P[2].m[1][1], (int)P[2].m[2][1],(int)P[3].m[1][1], (int)P[3].m[2][1]);
		}	
    	}
	}
  for (i = 0 ; i < TRIANGLE+1 ; i++) {
    free_dmatrix(P[i].m,1,P[i].l,1,P[i].c) ;
  }
}



int main() {

  Display *d ;
  Window w ;
  XEvent e ;
  int s ;

  unsigned int r, g, b ;
  double radius, a, c, theta, rho, delta_theta, delta_rho ;

  dmatrix_t E, G, C, L ; 
    
  dmat_alloc(&E,4,1) ; /* The centre of projection for the camera */
  E.m[1][1] = Ex ;
  E.m[2][1] = Ey ;
  E.m[3][1] = Ez ;
  E.m[4][1] = 1.0 ;
    
  dmat_alloc(&G,4,1) ; /* Point gazed at by camera */
  G.m[1][1] = Gx ;
  G.m[2][1] = Gy ;
  G.m[3][1] = Gz ;
  G.m[4][1] = 1.0 ;

  dmat_alloc(&L, 4, 1); /* location of light source */
  L.m[1][1] = Ex+20.0;
  L.m[2][1] = Ey+20.0;
  L.m[3][1] = Ez;
  L.m[4][1] = 1.0;

  dmat_alloc(&C,4,4) ;
  C = *build_camera_matrix(&E,&G) ; /* The camera matrix */

  delta_theta = 2.0*M_PI/48.0 ; 
  delta_rho = M_PI/18.0 ; 
  radius = 10.0 ;
  a = 5.0 ;
  c = 25.0 ; 

  r = g = b = 0 ;
  d = InitX(d,&w,&s) ;
  SetCurrentColorX(d,&(DefaultGC(d,s)),r,g,b) ;
  XNextEvent(d,&e) ;
  printf("Hello main\n");
  while (1) {
    XNextEvent(d,&e) ;
    if (e.type == Expose) {
      shade_wiremesh_sphere(d,w,s,&C, &L, radius,0.0,2.0*M_PI,0.0,M_PI,delta_theta,delta_rho) ;
      shade_wiremesh_torus(d,w,s,&C, &L, a,c,0.0,2.0*M_PI,0.0,2.0*M_PI,delta_theta,delta_theta) ;
      shade_wiremesh_cone(d,w,s,&C, &L, radius,0.0,2.0*M_PI,0.0,M_PI,delta_theta,delta_rho) ;
    }
    if(e.type == KeyPress)
      break ;
    if(e.type == ClientMessage)
      break ;
  }
  QuitX(d,w) ;
  free_dmatrix(E.m,1,E.l,1,E.c) ;
  free_dmatrix(G.m,1,G.l,1,G.c) ;
  free_dmatrix(C.m,1,C.l,1,C.c) ;
}
