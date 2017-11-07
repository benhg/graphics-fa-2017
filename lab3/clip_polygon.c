#include <FPT.h>   
#include "D2d_matrix.h"
#define WIDTH 600
#define HEIGHT WIDTH
int numpoints[10], numpolys[10];
double x[10][9000], y[10][9000];
double longest[10];
double s[10];
int psize[10][5000];
int con[10][5000][20];
double clipx[100], clipy[100];
double clipdx[100], clipdy[100];
int numclip;
double sectx, secty;
double red[10][5000], grn[10][5000], blu[10][5000];

void draw(int k) {
  int i, j;
  for(i = 0; i < numpolys[k]; i++) {
    G_rgb(red[k][i], grn[k][i], blu[k][i]);
    double polyx[9000], polyy[9000];
    for(j = 0; j < psize[k][i]; j++) {
      int point = con[k][i][j];
      polyx[j] = x[k][point];
      polyy[j] = y[k][point];
    }
    G_fill_polygon(polyx, polyy, psize[k][i]);
  }
}

int inside(int k, double xx, double yy, int c1, int c2) {
  int i;
  double avx, avy, avgres, res;
  avx = 0; avy = 0;
  for(i = 0; i < numclip; i++) {
    avx += clipx[i];
    avy += clipy[i];
  }
  avx = (avx / numclip);
  avy = (avy / numclip);

  avgres = (((clipy[c2] - clipy[c1]) * (avx - clipx[c1])) - ((clipx[c2] - clipx[c1]) * (avy - clipy[c1])));
  res = (((clipy[c2] - clipy[c1]) * (xx - clipx[c1])) - ((clipx[c2] - clipx[c1]) * (yy - clipy[c1])));
  if (avgres > 0 && res < 0) { return 1; }
  if (avgres < 0 && res > 0) { return 1; }
  else { return 0; }
}

void intersect(int k, double p1x, double p1y, double p2x, double p2y, int good1, int good2) {
  double a, b, c, p, q, r;
  a = (p2y - p1y);
  b = (p2x - p1x);
  c = ((b * p1y) - (a * p1x));
  p = (clipy[good2] - clipy[good1]);
  q = (clipx[good2] - clipx[good1]);
  r = ((q * clipy[good1]) - (p * clipx[good1]));
	
  sectx = (((r * b) - (c * q)) / ((a * q) - (p * b)));
  if(b == 0) {
    secty = (-1 * ((((-1 * p) * sectx) - r) / q));
  }  else {
    secty = (-1 * ((((-1 * a) * sectx) - c) / b));
  }
}

int clip_edge(int k, int size) {
  int l, a, b, z, good1, good2, newsize;
  for(l = 0; l < numclip; l++) {
    z = l + 1;  if (z == numclip) { z = 0; }
    double tempx[100], tempy[100];
    newsize = 0;

    for(a = 0; a < size; a++) {
      b = a + 1;  if (b == size) { b = 0; }
      good1 = inside(k, clipdx[a], clipdy[a], l, z);
      good2 = inside(k, clipdx[b], clipdy[b], l, z); 
      if(!good1 && !good2) {
        tempx[newsize] = clipdx[b];
        tempy[newsize] = clipdy[b];
        newsize++;
      } else if (!good1 && good2) {
        intersect(k, clipdx[a], clipdy[a], clipdx[b], clipdy[b], l, z);
        tempx[newsize] = sectx;
        tempy[newsize] = secty;
        newsize++;
      } else if (good1  && !good2) {
        intersect(k, clipdx[a], clipdy[a], clipdx[b], clipdy[b], l, z);
        tempx[newsize] = sectx;
        tempy[newsize] = secty;
        newsize++;

        tempx[newsize] = clipdx[b];
        tempy[newsize] = clipdy[b];
        newsize++;
      }
    }

    for(a = 0; a < newsize; a++) {
      clipdx[a] = tempx[a];
      clipdy[a] = tempy[a];
    }
    size = newsize;
  }
  return size;
}

void clip(int k) {
  int i, j, size, l;
  for(i = 0; i < numpolys[k]; i++) {
    for(j = 0; j < psize[k][i]; j++) {    
      clipdx[j] = x[k][con[k][i][j]];
      clipdy[j] = y[k][con[k][i][j]];
    }
    size = clip_edge(k, psize[k][i]);
    G_rgb(red[k][i], grn[k][i], blu[k][i]);
    G_fill_polygon(clipdx, clipdy, size);
  }

}

void stretch(double sx, double sy, int k) {
  double m[3][3], minv[3][3];
  D2d_make_identity(m); D2d_make_identity(minv);
  D2d_scale(m, minv,  sx, sy);
  D2d_mat_mult_points(x[k], y[k],  m, x[k], y[k], numpoints[k]);
}

void translate(double dx, double dy, int k) {
  double m[3][3], minv[3][3];
  D2d_make_identity(m); D2d_make_identity(minv);
  D2d_translate(m, minv, dx, dy);
  D2d_mat_mult_points(x[k], y[k],  m, x[k], y[k], numpoints[k]);
}

void rotate(int k, double rads) {
  translate(-300, -300, k);
  double m[3][3], minv[3][3];
  D2d_make_identity(m); D2d_make_identity(minv);
  D2d_rotate(m, minv,  rads);
  D2d_mat_mult_points(x[k], y[k],  m, x[k], y[k], numpoints[k]);
  translate(300, 300, k);
}

void draw_clip_window() {
  int i;
  for(i = 0; i < numclip - 1; i++) {
    G_line(clipx[i], clipy[i], clipx[i + 1], clipy[i + 1]);
  }
  G_line(clipx[numclip - 1], clipy[numclip - 1], clipx[0], clipy[0]);
}

void set_globals(FILE *q, int key, int argc,  char **argv){
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


int main(int argc,  char **argv) {
  FILE *q;
  int key;

  G_init_graphics(WIDTH,HEIGHT);

  G_rgb(0,0,0) ;
  G_clear() ;
  set_globals(q, key, argc, argv);
  
  draw(1);

  G_rgb(1,1,1) ;
  G_fill_rectangle(0,0,50,50);
  G_rgb(1,1,0) ;

  double p[2];
  numclip = 0;
  G_wait_click(p) ; 
  clipx[numclip] = p[0]; clipy[numclip] = p[1];
  G_point(p[0], p[1]);
  numclip++;
  while (1) {
    G_wait_click(p) ; 
    if ((p[0] < 50) && (p[1] < 50)) { break; }
    clipx[numclip] = p[0]; clipy[numclip] = p[1];
    G_line(clipx[numclip - 1], clipy[numclip - 1], clipx[numclip], clipy[numclip]);
    numclip++;
  }
  G_line(clipx[numclip - 1], clipy[numclip - 1], clipx[0], clipy[0]);
 
  int lastKey = 0;
  while (1) {
    key = G_wait_key();
    G_rgb(0, 0, 0);
    G_clear();
    if (key == 'x') {break;}
    if (key == lastKey) {
      rotate(key - 48, .05);
      clip(key - 48);
      draw_clip_window();
    } else {
      clip(key - 48);
      draw_clip_window();
    }
    lastKey = key;
  }
}




