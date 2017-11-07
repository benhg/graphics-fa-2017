#include <FPT.h>   
#include "D3d_matrix.h"

int numpoints[10], numpolys[10];
int flipflag;
double x[10][9000], y[10][9000], z[10][9000];
double xbb[10][9000], ybb[10][9000], zbb[10][9000];
double avgx[10], avgy[10], avgz[10];
double longest[10];
double s[10];
int psize[10][5000];
int con[10][5000][20];
double halfangle;

void printxy(int k) {
  int i;
  for (i = 0; i < numpoints[k]; i++) {
    printf("x: %f y: %f z: %f \n", x[k][i], y[k][i], z[k][i]);
  }
  printf("\n");
}

void printbb(int k) {
  int i;
  for (i = 0; i < numpoints[k]; i++) {
    printf("x: %f y: %f \n", xbb[k][i], ybb[k][i]);
  }
  printf("\n");
}

void prep(int k) {
  int i, j;
  double h;
  h = tan(halfangle);
  for(i = 0; i < numpoints[k]; i++) {
      xbb[k][i] = ((x[k][i] / z[k][i]) * (300 / h) + 300);
      ybb[k][i] = ((y[k][i] / z[k][i]) * (300 / h) + 300);
  }
}

void draw(int k) {
  prep(k);
  int i, j;
  for(i = 0; i < numpolys[k]; i++) {
    double polyx[9000], polyy[9000], n[3], v[3], a[3], b[3], dotprod;
    int point, point1, point2, point3;
    point1 = con[k][i][0];
    point2 = con[k][i][1];
    point3 = con[k][i][(psize[k][i] - 1)];
    a[0] = (x[k][point2] - x[k][point1]); a[1] = (y[k][point2] - y[k][point1]); a[2] = (z[k][point2] - z[k][point1]);
    b[0] = (x[k][point3] - x[k][point1]); b[1] = (y[k][point3] - y[k][point1]); b[2] = (z[k][point3] - z[k][point1]);
    v[0] = x[k][point1]; v[1] = y[k][point1]; v[2] = z[k][point1];
    D3d_x_product(n, a, b);
    dotprod = ((v[0] * n[0]) + (v[1] * n[1]) + (v[2] * n[2]))*flipflag;
    if (dotprod > 0) {
      for(j = 0; j < psize[k][i]; j++) {
        point = con[k][i][j];
        polyx[j] = xbb[k][point];
        polyy[j] = ybb[k][point];
     	}
    G_polygon(polyx, polyy, psize[k][i]);
    }
  }
}

void stretch(double sx, double sy, double sz, int k) {
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_scale(m, minv,  sx, sy, sz);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
}

void translate(double dx, double dy, double dz, int k) {
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_translate(m, minv, dx, dy, dz);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
}

void rotx(int k, double rads) {
  int i;
  double avgy, avgz;
  avgy = 0; avgz = 0;
  for(i = 0; i < numpoints[k]; i++) {
    avgy += y[k][i];
    avgz += z[k][i];
  }
  avgy = (avgy / numpoints[k]);
  avgz = (avgz / numpoints[k]);
  translate(0, -avgy, -avgz, k);//move z?
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_x(m, minv,  rads);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
  translate(0, avgy, avgz, k);
}

void roty(int k, double rads) {
  int i;
  double avgx, avgz;
  avgx = 0; avgz = 0;
  for(i = 0; i < numpoints[k]; i++) {
    avgx += x[k][i];
    avgz += z[k][i];
  }
  avgx = (avgx / numpoints[k]);
  avgz = (avgz / numpoints[k]);
  translate(-avgx, 0, -avgz, k);
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_y(m, minv,  rads);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
  translate(avgx, 0, avgz, k);
}

void rotz(int k, double rads) {
  int i;
  double avgx, avgy;
  avgx = 0; avgy = 0;
  for(i = 0; i < numpoints[k]; i++) {
    avgx += x[k][i];
    avgy += y[k][i];
  }
  avgx = (avgx / numpoints[k]);
  avgy = (avgy / numpoints[k]);
  translate(-avgx, -avgy, 0, k);
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_z(m, minv,  rads);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
  translate(avgx, avgy, 0, k);
}

int main(int argc,  char **argv) {
  FILE *q;
  int i, j, k, key;
  
  halfangle = .5;

  G_init_graphics(600,600);
  G_rgb(0,0,0);
  G_clear();

  for(k = 1; k < argc; k++) {
    q = fopen(argv[k], "r");
    if (q == NULL) {
      printf("can't open file \n");
      exit(0);
    }

    fscanf(q, "%d", &numpoints[k]);
 
    for(i = 0; i < numpoints[k]; i++) {
      fscanf(q, "%lf %lf %lf", &x[k][i], &y[k][i], &z[k][i]); 
    }

    fscanf(q, "%d", &numpolys[k]);

    for(i = 0; i < numpolys[k]; i++) {
      fscanf(q, "%d", &psize[k][i]);
      for(j = 0; j < psize[k][i]; j++) {
        fscanf(q, "%d", &con[k][i][j]);
      }
    }
  }
  
  G_rgb(1, 0, 0);
  draw(1);

  G_wait_key();
 
  k = 1;
  flipflag = 1;
  while (0 == 0) {
    key = G_wait_key();
    G_rgb(0, 0, 0);
    G_clear(k);
    G_rgb(1, 0, 1);
    if (key == 'x') {
      rotx(k, .1);
      draw(k);
    } else if (key == 'c') {
      roty(k, .1);
      draw(k);
    } else if (key == 'z') {
      rotz(k, .1);
      draw(k);
    }else if (key == 'X') {
      rotx(k, -.1);
      draw(k);
    } else if (key == 'C') {
      roty(k, -.1);
      draw(k);
    } else if (key == 'Z') {
      rotz(k, -.1);
      draw(k);
    }
    else if (key == 'd' || key == 'D') {
      translate(1, 0, 0, k);
      draw(k);
    } else if (key == 'a' || key == 'A') {
      translate(-1, 0, 0, k);
      draw(k);
    } else if (key == 's' || key == 'S') {
      translate(0, -1, 0, k);
      draw(k);
    } else if (key == 'w' || key == 'W') {
      translate(0, 1, 0, k);
      prep(k);
      draw(k);
    } else if (key == 'q' || key == 'Q') {
      translate(0, 0, 1, k);
      draw(k);
    } else if (key == 'e' || key == 'E') {
      translate(0, 0, -1, k);
      draw(k);
    }else if(key=='f' || key == 'F'){
      flipflag *= -1;
      draw(k);
    }else if(key=='y'){
      break;
    }
    else if (key > 47 && key < 58) {
      draw(k);
      k = (key - 48);
    }
  }
  
}
