#include <FPT.h>   
#include "D3d_matrix.h"
#define TRUE 1
#define FALSE !TRUE
xs
double x[100], y[100], z[100];


void translate(double dx, double dy, double dz, int numpts) {
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_translate(m, minv, dx, dy, dz);
  D3d_mat_mult_points(x, y, z, m, x, y, z, numpts);
}

void rotate(double rads, int numpts) {
  translate(-300, 0, 0, numpts);
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_y(m, minv,  rads);
  D3d_mat_mult_points(x, y, z, m, x, y, z, numpts);
  translate(300, 0, 0, numpts);
}

int main(int argc,  char **argv) {
  FILE *f;
  f = fopen(argv[1], "w");
  if (f == NULL) {
    printf("can't open file \n");
    exit(0);
  }

  int numpts, n, i, j;
  double P[2], rads;
  printf("Input n of passes\n");
  scanf("%d", &n);

  G_init_graphics(600, 600);
  G_rgb(0,0,0);
  G_clear();
  G_rgb(1,0,1);
  G_fill_rectangle(0, 0, 50, 50);
  
  G_rgb(1,0,0);
  G_line(300, 0, 300, 600);
  G_wait_click(P);
  x[0] = P[0]; y[0] = P[1]; z[0] = 0;
  numpts = 0;
  G_fill_circle(P[0], P[1], 2);
  
  while(TRUE) {
    G_wait_click(P);
    if(P[0] < 50 && P[1] < 50) { break; }
    G_fill_circle(P[0], P[1], 2);
    G_line(P[0], P[1], x[numpts], y[numpts]);
    numpts++;
    x[numpts] = P[0]; y[numpts] = P[1]; z[0] = 0;
  }
  numpts++;
  
  rads = ((2 * M_PI) / n);
  fprintf(f, "%d\n", (n * numpts));
  for(i = 1; i <= n; i++) {
    rotate(rads, numpts);
    for(j = 0; j < numpts; j++) {
      fprintf(f, "%f %f %f\n", x[j], y[j], z[j]);
    }
  }
  
  fprintf(f, "%d\n", n * (numpts - 1));
  for(i = 0; i < n; i++) {
    for(j = 0; j < (numpts - 1); j++) {
      fprintf(f, "4 ");
      fprintf(f, "%d ", ((i * numpts) + j));
      fprintf(f, "%d ", ((i * numpts) + j + 1));
      if(((i * numpts) + j + numpts + 1) > (n * numpts)) {
        fprintf(f, "%d ", j + 1);
        fprintf(f, "%d\n", j);
      } else {
        fprintf(f, "%d ", ((i * numpts) + j + numpts + 1));
        fprintf(f, "%d\n", ((i * numpts) + j + numpts));
      }
    }
  }
  
}

















