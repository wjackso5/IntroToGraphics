#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define M_PI (3.14159)

void Bresenham(int, int, int, int, Display*, Window, GC);
	
int main(){
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
		//draw the assignment sample
		double dt = 2.0*M_PI/200.0 ;
		for (double t = 0.0 ; t < 2.0*M_PI;) {
	  		int x1 = 256 + (int)100.0*(1.5*cos(t) - cos(13.0*t)) ;
	  		int y1 = 256 + (int)100.0*(1.5*sin(t) - sin(13.0*t)) ;
	  		t += dt ;
	  		int x2 = 256 + (int)100.0*(1.5*cos(t) - cos(13.0*t)) ;
	  		int y2 = 256 + (int)100.0*(1.5*sin(t) - sin(13.0*t)) ;
	  		Bresenham(x1,y1,x2,y2,disp,window,gc); //This function draws points
		}
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
	return 0;
}


void Bresenham(int x1, int y1, int x2, int y2, Display *disp, Window window, GC gc){
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
                     
