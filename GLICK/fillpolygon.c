#include <FPT.h>
#define TRUE 1
#define FALSE !TRUE

double xpt[100];
int GRAPHICS_HEIGHT;
int GRAPHICS_WIDTH;
GRAPHICS_WIDTH = 1000,
    GRAPHICS_HEIGHT = 1000;

void selsort(double * xpoi, int poi) {
    int i, j;
    double swap;
    for (i = 0; i < poi - 1; i++) {
        for (j = 0; j < poi - i - 1; j++) {
            if (xpoi[j] > xpoi[j + 1]) {
                swap = xpoi[j];
                xpoi[j] = xpoi[j + 1];
                xpoi[j + 1] = swap;
            }
        }
    }
}

void fill_polygon(double * p_x, double * p_y, int n) {
    int i, j, poi, k;
    double y;
    for (y = 0; y < GRAPHICS_HEIGHT; y++) {
        poi = 0;
        k = n - 1;
        for (j = 0; j < n; j++) {
            if ((p_y[j] < y && p_y[k] >= y) || (p_y[k] < y && p_y[j] >= y)) {
                xpt[poi++] = ((p_x[k] - p_x[j]) * (y - p_y[k])) / (p_y[k] - p_y[j]) + p_x[k];
            }
            k = j;
        }
        selsort(xpt, poi);
        for (i = 0; i < poi; i += 2) {
	  G_line(xpt[i], y, xpt[i + 1], y);
        }
    }
}

void draw_rectangles() {
    G_fill_rectangle(0, 0, 50, 50);
    G_fill_rectangle(GRAPHICS_WIDTH - 50, 0, GRAPHICS_WIDTH, 50);
    G_fill_rectangle(GRAPHICS_WIDTH - 50, GRAPHICS_HEIGHT - 50, GRAPHICS_WIDTH, GRAPHICS_HEIGHT);
    G_fill_rectangle(0, GRAPHICS_HEIGHT - 50, 50, GRAPHICS_HEIGHT);
}
int detect(double *p){
  return 1;
}

int main() {
    G_init_graphics(GRAPHICS_WIDTH, GRAPHICS_HEIGHT);
    G_rgb(0, 0, 0);
    draw_rectangles();
    int i, n, j, k, poi, FLAG;
    double y;
    double p[2];
    double p_x[100];
    double p_y[100];
    i = 0;
    n = 0;
    
    while (TRUE) {

        G_wait_click(p);
        if (p[0] > GRAPHICS_WIDTH - 50 && p[1] > GRAPHICS_HEIGHT - 50) {
            G_rgb(0, 0, 255);
            //fill_polygon(p_x, p_y, n);
	    if(detect(p)){
	      G_rgb(0,255,0);
	    }else{
	      G_rgb(255,0,0);
	    }
	    G_fill_rectangle(p[0] - 2, p[1] - 2, 4, 4);
	    FLAG = 1;
            G_rgb(0, 0, 0);
        }
        if (p[0] < 50 && p[1] < 50) {
            break;
        }
        if ((p[0] < 50) && (p[1] > GRAPHICS_HEIGHT - 50)) {
            G_rgb(255, 0, 0);
            G_polygon(p_x, p_y, n);
	    FLAG = 1;
            G_rgb(0, 0, 0);
        }
        if (!(p[0] < GRAPHICS_WIDTH - 50) && (p[1] < 50)) {
            G_rgb(255, 255, 255);
            G_clear();
            G_rgb(0, 0, 0);
            i = 0;
            n = 0;
            draw_rectangles();
        }
        if (!(p[0] > GRAPHICS_WIDTH - 50) && !(p[1] > GRAPHICS_HEIGHT - 50) && !(p[0] > GRAPHICS_WIDTH - 50) && !(p[1] < 50)) {
	  if (FLAG==1){
	     i = 0;
            n = 0;
	    FLAG = 0;
	  }
	    G_fill_rectangle(p[0] - 2, p[1] - 2, 4, 4);
            p_x[i] = p[0];
            p_y[i] = p[1];
            i++;
            n++;
        }
    }

}
