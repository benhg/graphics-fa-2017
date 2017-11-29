#include <FPT.h>   
#include <D3d_matrix.h>

int numpoints[10], numpolys[10], flag, numobs, totpolys;
double x[10][9000], y[10][9000], z[10][9000];
double xbb[10][9000], ybb[10][9000], zbb[10][9000];
double avgx[10], avgy[10], avgz[10];
double longest[10];
double s[10];
int psize[10][5000];
int con[10][5000][20];
double halfangle;
typedef struct {int objnum; int polynum; double dist;} THING;
THING things[50000];

int compare (const void *p, const void *q) {
  THING *a, *b;
  a = (THING*)p;
  b = (THING*)q;
  if  (((*a).dist) < ((*b).dist)) return 1;
  else if (((*a).dist) > ((*b).dist)) return -1;
  else return 0;
}

void sortit() {
  int i, j, q, point;
  double d;
  q = 0;
  for(i = 1; i <= numobs; i++) {
    for(j = 0 ; j < numpolys[i]; j++) {
      point = con[i][j][0];
      d = sqrt((pow(x[i][point], 2)) + (pow(y[i][point], 2)) + (pow(z[i][point], 2)));
      things[q].objnum = i; things[q].polynum = j; things[q].dist = d;
      q++;
    }
  }  
  qsort(things, totpolys, sizeof(THING), compare);
}

void prep() {
  int i, j;
  double h;
  h = tan(halfangle);
  for(int j = 1 ; j <= numobs; j++){
    for(i = 0; i < numpoints[j]; i++) {
      xbb[j][i] = ((x[j][i] / z[j][i]) * (300 / h) + 300);
      ybb[j][i] = ((y[j][i] / z[j][i]) * (300 / h) + 300);
    }
  }
}

void normalize(double a[]) {
  double dist;
  dist = sqrt((pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2)));
  a[0] = (a[0] / dist); a[1] = (a[1] / dist); a[2] = (a[2] / dist);
}

double color(double intensity, double actrgb[], double ambient, double diffusemax) {
  double rgb[3], inhrgb; rgb[0] = 1; rgb[1] = 0; rgb[2] = 0;
  inhrgb = ambient + diffusemax;
  if (intensity < inhrgb) {
    inhrgb = (intensity / inhrgb);
    actrgb[0] = inhrgb * rgb[0]; actrgb[1] = inhrgb * rgb[1]; actrgb[2] = inhrgb * rgb[2];
  } else if (intensity > inhrgb) {
    actrgb[0] = ((1 - rgb[0]) * ((1 - intensity) / (1 - inhrgb))) + rgb[0];
    actrgb[1] = ((1 - rgb[1]) * ((1 - intensity) / (1 - inhrgb))) + rgb[1];
    actrgb[2] = ((1 - rgb[2]) * ((1 - intensity) / (1 - inhrgb))) + rgb[2];
  } else {
    actrgb[0] = inhrgb * rgb[0]; actrgb[1] = inhrgb * rgb[1]; actrgb[2] = inhrgb * rgb[2];
  }
}

void draw() {
	int i, j, point, point1, point2, point3, k;
	sortit();
	prep();
	double light[3]; light[0] = 100; light[1] = 200; light[2] = 0;
	double diffusemax, ambient, intensity, specpow, specular, diffuse;
	diffusemax = .6; ambient = .2; specpow = 30;
  for(i = 0; i < totpolys; i++) {  
    double polyx[9000], polyy[9000], n[3], a[3], b[3], l[3], e[3], r[3], nde, ndl, dist, actrgb[3];
    for(j = 0; j < psize[things[i].objnum][things[i].polynum]; j++) {
      point = con[things[i].objnum][things[i].polynum][j];
      polyx[j] = xbb[things[i].objnum][point];
      polyy[j] = ybb[things[i].objnum][point];
    }    
    
    k = (things[i].objnum);
    point1 = con[k][things[i].polynum][0];
    point2 = con[k][things[i].polynum][1];
    point3 = con[k][things[i].polynum][2];
    a[0] = (x[k][point2] - x[k][point1]); a[1] = (y[k][point2] - y[k][point1]); a[2] = (z[k][point2] - z[k][point1]);
    b[0] = (x[k][point3] - x[k][point2]); b[1] = (y[k][point3] - y[k][point2]); b[2] = (z[k][point3] - z[k][point2]);    
    l[0] = (light[0] - x[k][point1]); l[1] = (light[1] - y[k][point1]); l[2] = (light[2] - z[k][point1]); normalize(l);
    e[0] = (0 - x[k][point1]); e[1] = (0 - y[k][point1]); e[2] = (0 - z[k][point1]); normalize(e);
    D3d_x_product(n, a, b); normalize(n);    
    ndl = ((n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2]));
    nde = ((n[0] * e[0]) + (n[1] * e[1]) + (n[2] * e[2]));
    if(nde < 0 && ndl < 0) { n[0] = (-1 * n[0]); n[1] = (-1 * n[1]); n[2] = (-1 * n[2]); }
    ndl = ((n[0] * l[0]) + (n[1] * l[1]) + (n[2] * l[2]));
    nde = ((n[0] * e[0]) + (n[1] * e[1]) + (n[2] * e[2]));    
    if((nde * ndl) < 0) {
      intensity = ambient; 
    } else {
      r[0] = (((2 * ndl) * n[0]) - l[0]);
      r[1] = (((2 * ndl) * n[1]) - l[1]);
      r[2] = (((2 * ndl) * n[2]) - l[2]); normalize(r);    
      diffuse = (diffusemax * ndl);
      double edotr = (e[0] * r[0]) + (e[1] * r[1]) + (e[2] * r[2]) ;
      if (edotr < 0) { edotr = 0 ; }      
      specular = ((pow(edotr, specpow)) * (1 - ambient - diffusemax));    
      intensity = ambient + diffuse + specular;
    }
    
    color(intensity, actrgb, ambient, diffusemax);    
    G_rgb(actrgb[0], actrgb[1], actrgb[2]);
    //G_rgb(intensity, intensity, intensity);
    G_fill_polygon(polyx, polyy, psize[things[i].objnum][things[i].polynum]);
    //G_rgb(0,0,0);
    //G_polygon(polyx, polyy, psize[things[i].objnum][things[i].polynum]);
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
  translate(0, -avgy, -avgz, k);
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
  flag = 1;
  numobs = argc - 1;

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
      fscanf(q, "%lf %lf %lf", &x[k][i], &y[k][i], &z[k][i]); //writes to the array
    }

    fscanf(q, "%d", &numpolys[k]);

    for(i = 0; i < numpolys[k]; i++) {
      totpolys++;
      fscanf(q, "%d", &psize[k][i]);
      for(j = 0; j < psize[k][i]; j++) {
        fscanf(q, "%d", &con[k][i][j]);
      }
    }
  }
  
  
  //special translation for writefile
  double avgx, avgy, avgz;
  for(i = 0; i < numpoints[1]; i++) {
    avgx += x[1][i];
    avgy += y[1][i];
    avgz += z[1][i];
  }
  avgy = -(avgy / numpoints[1]);
  avgx = -(avgx / numpoints[1]);
  avgz = 500 + (avgz / numpoints[1]);
  translate(avgx, avgy, avgz, 1);
  
  draw();
  G_wait_key(); 
  k = 1;
  while (0 == 0) {
    key = G_wait_key();
    G_rgb(0, 0, 0);
    G_clear();
    G_rgb(1, 0, 0);
     if (key == 'x') {
      rotx(k, .1);
      draw();
    } else if (key == 'c') {
      roty(k, .1);
      draw();
    } else if (key == 'z') {
      rotz(k, .1);
      draw();
    }else if (key == 'X') {
      rotx(k, -.1);
      draw();
    } else if (key == 'C') {
      roty(k, -.1);
      draw();
    } else if (key == 'Z') {
      rotz(k, -.1);
      draw();
    }
    else if (key == 'd' || key == 'D') {
      translate(1, 0, 0, k);
      draw();
    } else if (key == 'a' || key == 'A') {
      translate(-1, 0, 0, k);
      draw();
    } else if (key == 's' || key == 'S') {
      translate(0, -1, 0, k);
      draw();
    } else if (key == 'w' || key == 'W') {
      translate(0, 1, 0, k);
      prep();
      draw();
    } else if (key == 'q' || key == 'Q') {
      translate(0, 0, 1, k);
      draw();
    } else if (key == 'e' || key == 'E') {
      translate(0, 0, -1, k);
      draw();
    }else if(key=='f' || key == 'F'){
      draw();
    }else if(key=='y'){
      break;
    }
    else if (key > 47 && key < 58) {
      draw();
      k = (key - 48);
    }
  }
  
}







