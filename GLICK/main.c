#include <FPT.h>   
#include <D2d_matrix.h>
#define WIDTH 600
#define HEIGHT 600
#define TRUE 1
#define FALSE !TRUE

int numpoints[10], numpolys[10];
double x[10][9000], y[10][9000];
double s[10];
int psize[10][5000];
int con[10][5000][20];
double red[10][5000], grn[10][5000], blu[10][5000];

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


void stretch(double sx, double sy, int k) {//what is th inv for??
  double m[3][3], minv[3][3];
  D2d_make_identity(m);
  D2d_make_identity(minv);
  D2d_scale(m, minv,  sx, sy);
  D2d_mat_mult_points(x[k], y[k],  m, x[k], y[k], numpoints[k]);
}

void translate(double dx, double dy, int k) {
  double m[3][3], minv[3][3];
  D2d_make_identity(m);
  D2d_make_identity(minv);
  D2d_translate(m, minv, dx, dy);
  D2d_mat_mult_points(x[k], y[k],  m, x[k], y[k], numpoints[k]);
}

void rotate(double radians, int k) {
  translate(-(HEIGHT/2), -(HEIGHT/2), k);
  double m[3][3], minv[3][3];
  D2d_make_identity(m); D2d_make_identity(minv);
  D2d_rotate(m, minv,  radians);
  D2d_mat_mult_points(x[k], y[k],  m, x[k], y[k], numpoints[k]);
  translate((HEIGHT/2), (HEIGHT/2), k);
}

void load_points(FILE *q, int key, int argc,  char **argv){
  int i, j, k;
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
	printf("%lf\n", distance);
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
    int key;
    double radians = 0.05;
    G_init_graphics(HEIGHT, WIDTH);

    G_rgb(0, 0, 0);
    G_clear();
    load_points(q, key,argc, argv);
    while (TRUE) {
        key = G_wait_key();
        key = key - 48;
        G_rgb(0, 0, 0);
        G_clear();
        if (key == 72) {
            break;
        }
        draw(key);
        rotate(radians, key);
    }

}
