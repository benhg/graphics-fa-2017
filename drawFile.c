#include <FPT.h>
#define TRUE 1
#define FALSE !TRUE
#define HEIGHT 800
#define WIDTH 800

int numpoints[10], numpolys[10];
double x[10][9000], y[10][9000];
int psize[10][5000]; //size of polyogns
int con[10][5000][20]; //rc, row column
double red[10][5000], grn[10][5000], blu[10][5000];

void translate(double dx, double dy, int onum) {
    int i;
    for (i = 0; i < numpoints[onum]; i++) {
        x[onum][i] += dx;
        y[onum][i] += dy;
    }
}

void stretch(double sx, double sy, int onum) {
    int i;
    for (i = 0; i < numpoints[onum]; i++) {
        x[onum][i] *= sx;
        y[onum][i] *= sy;
    }
}

void rotate(double theta, int onum) {
    int i;
    double x_n, y_n;
    translate(-(WIDTH/2),-(HEIGHT/2),onum);
    for (i = 0; i < numpoints[onum]; i++) {
        x_n = x[onum][i];
        y_n = y[onum][i];
        x[onum][i] = (x_n * cos(theta)) - (y_n * sin(theta));
        y[onum][i] = (y_n * cos(theta)) + (x_n * sin(theta));
    }
    translate((WIDTH/2),(HEIGHT/2),onum);
}

void draw(int k) {
    int i, j;
    for (i = 0; i < numpolys[k]; i++) {
        G_rgb(red[k][i], grn[k][i], blu[k][i]);
        double polyx[9000], polyy[9000];
        for (j = 0; j < psize[k][i]; j++) {
            int point = con[k][i][j];
            polyx[j] = x[k][point];
            polyy[j] = y[k][point];
        }
        G_fill_polygon(polyx, polyy, psize[k][i]);
    }
}

void load_points(FILE * q, int i,int j,int k,int key, int argc,  char ** argv){
    for (k = 1; k < argc; k++) {
        q = fopen(argv[k], "r");
        if (q == NULL) {
            printf("can't open file \n");
            exit(0);
        }

        fscanf(q, "%d", & numpoints[k]);

        for (i = 0; i < numpoints[k]; i++) {
            fscanf(q, "%lf %lf", & x[k][i], & y[k][i]); //writes to the array
        }
	double totalx, totaly, meanx, meany, distance, first, second, diff;

	distance = 0;
	for (i = 0; i < numpoints[k]; i++) {
          for(j=0;j<numpoints[k];j++){
	    first = x[k][i];
	    second = x[k][j];
	    diff = second-first;
	    if(diff>=distance){
	      distance = diff;
	    }
	  }
        }
	printf("%lf/n", distance);
	
       	stretch(((WIDTH/(2*distance))),((WIDTH/(2*distance))), k);
	
	totalx = 0;
	totaly = 0;
	for (i = 0; i < numpoints[k]; i++) {
          totalx += x[k][i];
	  totaly += y[k][i];
        }
	meanx = (totalx)/(numpoints[k])-(WIDTH/2);
	meany = (totaly)/(numpoints[k])-(HEIGHT/2);
	translate(-meanx, -meany, k);

        fscanf(q, "%d", & numpolys[k]);

        for (i = 0; i < numpolys[k]; i++) {
            fscanf(q, "%d", & psize[k][i]);
            for (j = 0; j < psize[k][i]; j++) {
                fscanf(q, "%d", & con[k][i][j]);
            }
        }

        for (i = 0; i < numpolys[k]; i++) {
            fscanf(q, "%lf %lf %lf", & red[k][i], & grn[k][i], & blu[k][i]);
        }
    }
}

int main(int argc, char ** argv) {
    FILE * q;
    int i, j, k, key;
    G_init_graphics(HEIGHT, WIDTH);

    G_rgb(0, 0, 0);
    G_clear();
    load_points(q,i,j,k,key,argc, argv);
    while (1) {
        key = G_wait_key();
        key = key - 48;
        G_rgb(0, 0, 0);
        G_clear();
        if (key == 72) {
            break;
        }
        draw(key);
        rotate(0.2, key);
    }

}
