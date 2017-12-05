#include <FPT.h>   
#include "D3d_matrix.h"

int numpoints[10], numpolys[10], flag, numobs, totpolys;
double x[10][9000], y[10][9000], z[10][9000];
double xbb[10][9000], ybb[10][9000], zbb[10][9000];
double avgx[10], avgy[10], avgz[10];
double longest[10];
double s[10];
double eye[3], coi[3], up[3];
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

int planeclip(double a, double b, double c, double d, double polyx[], double polyy[], double polyz[],
              int numpts, double resx[], double resy[], double resz[]) {
  int num, i, j;
  double x1, y1, z1, x2, y2, z2, x21, y21, z21, den, t, xintsct, yintsct, zintsct, s1, s2;
  
  num = 0;
  for (i = 0; i < numpts; i++) {
    j = (i + 1) % numpts;
    x1 = polyx[i]; y1 = polyy[i]; z1 = polyz[i];
    x2 = polyx[j]; y2 = polyy[j]; z2 = polyz[j];
    s1 = (a*x1 + b*y1 + c*z1 + d);
    s2 = (a*x2 + b*y2 + c*z2 + d);
    if ((s1 >= 0) && (s2 >= 0)) { 
      //out do nothing    
    } else if ((s1 < 0) && (s2 < 0)) {
      resx[num] = x2; resy[num] = y2; resz[num] = z2; num++;
    } else {      
      // one is in, the other out, so find the intersection
      x21 = x2 - x1; y21 = y2 - y1; z21 = z2 - z1;
      den = a*x21 + b*y21 + c*z21;
      if (den == 0) continue; 
      t = -(a*x1 + b*y1 + c*z1 + d) / den;
      xintsct = x1 + t*x21;
      yintsct = y1 + t*y21;
      zintsct = z1 + t*z21;
      if (s1 < 0) { 
        // in to out
        resx[num] = xintsct; resy[num] = yintsct; resz[num] = zintsct; num++;
      } else  {
        // out to in
        resx[num] = xintsct; resy[num] = yintsct; resz[num] = zintsct; num++;
        resx[num] = x2     ; resy[num] = y2     ; resz[num] = z2     ; num++;
      }
    }
  }
  return num;
}

int clip(double polyx[], double polyy[], double polyz[], int numpts) {
  double hith = 5; double yon = 100; double cent = ((hith + yon) / 2);
  double a, b, c, d, yh, hh, wx[6], wy[6], wz[6], nx[9000], ny[9000], nz[9000], cross1[3], cross2[3], n[3];
  int k, m, i;
  
  yh = (tan(halfangle) * yon);
  hh = (tan(halfangle) * hith);
  wx[0] = yh; wx[1] = yh; wx[2] = -yh; wx[3] = -yh;
  wy[0] = yh; wy[1] = -yh; wy[2] = -yh; wy[3] = yh;
  wz[0] = yon; wz[1] = yon; wz[2] = yon; wz[3] = yon;
    
  for(k = 0 ; k < 4; k++) { 
    m = k+1; if(m == 4) { m = 0; }
    cross1[0] = wx[k]; cross1[1] = wy[k]; cross1[2] = wz[k];
    cross2[0] = wx[m]; cross2[1] = wy[m]; cross2[2] = wz[m];
    D3d_x_product(n, cross1, cross2);
    a = n[0]; b = n[1]; c = n[2]; d = 0;
    if (a*0 + b*0 + c*cent + d > 0) {
     a = -a; b = -b; c = -c; d = -d;
    }
  
    numpts = planeclip(a, b, c, d, polyx, polyy, polyz, numpts, nx, ny, nz);    
    for (i = 0 ; i < numpts ; i++) {
      polyx[i] = nx[i];   polyy[i] = ny[i]; polyz[i] = nz[i];  
    }
  }
  //front and back clips
  a = 0; b = 0; c = 1; d = -yon;
  if (a*0 + b*0 + c*cent + d > 0) {
   a = -a; b = -b; c = -c; d = -d;;
  }
  numpts = planeclip(a, b, c, d, polyx, polyy, polyz, numpts, nx, ny, nz);    
  for (i = 0 ; i < numpts ; i++) {
    polyx[i] = nx[i];   polyy[i] = ny[i]; polyz[i] = nz[i];  
  }
  //again for the front
  a = 0; b = 0; c = 1; d = -hith;
  if (a*0 + b*0 + c*cent + d > 0) {
   a = -a; b = -b; c = -c; d = -d;
  }
  numpts = planeclip(a, b, c, d, polyx, polyy, polyz, numpts, nx, ny, nz);    
  for (i = 0 ; i < numpts ; i++) {
    polyx[i] = nx[i];   polyy[i] = ny[i]; polyz[i] = nz[i];  
  }
  return numpts;
}

void draw() {
	int i, j, point, point1, point2, point3, k, nsize, numpts;
	sortit();
	double h = tan(halfangle);
	double light[3]; light[0] = 100; light[1] = 200; light[2] = 0;
	double diffusemax, ambient, intensity, specpow, specular, diffuse; diffusemax = .5; ambient = .2; specpow = 30;
  for(i = 0; i < totpolys; i++) {  
    double polyx[9000], polyy[9000], polyz[9000], n[3], a[3], b[3], l[3], e[3], r[3], nde, ndl, dist, actrgb[3];
    
    for(j = 0; j < psize[things[i].objnum][things[i].polynum]; j++) {
      point = con[things[i].objnum][things[i].polynum][j];      
      polyx[j] = x[things[i].objnum][point];
      polyy[j] = y[things[i].objnum][point];
      polyz[j] = z[things[i].objnum][point];
    }
    numpts = clip(polyx, polyy, polyz, psize[things[i].objnum][things[i].polynum]);    
    for(j = 0; j < numpts; j++) {
      polyx[j] = ((polyx[j] / polyz[j]) * (300 / h) + 300);
      polyy[j] = ((polyy[j] / polyz[j]) * (300 / h) + 300);      
    }
    
    //color
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
    G_fill_polygon(polyx, polyy, numpts);
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

void eye_move(double dx, double dy, double dz){
     double m[4][4], minv[4][4];
     int k;
     for(k=0;k<numobs+1;k++){
       translate(dx,dy,dz,k);
     }
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

void eye_rotx(int k, double rads) {
  int i;
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_x(m, minv,  rads);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
}

void eye_roty(int k, double rads) {
  int i;
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_y(m, minv,  rads);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
}

void eye_rotz(int k, double rads) {
  int i;
  double m[4][4], minv[4][4];
  D3d_make_identity(m); D3d_make_identity(minv);
  D3d_rotate_z(m, minv,  rads);
  D3d_mat_mult_points(x[k], y[k], z[k], m, x[k], y[k], z[k], numpoints[k]);
}



void eye_rotate(double radsx, double radsy, double radsz){
     int k;
     for(k=0;k<numobs+1;k++){
       eye_rotx(k, radsx);
       eye_roty(k, radsy);
       eye_rotz(k, radsz);
     }

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
  //translate(avgx, avgy, avgz, 1);
  
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
    }else if (key == 'i') {
      eye_move(0,1,0);
      draw();
    } else if (key == 'j') {
       eye_move(-1,0,0);
       draw();
    } else if (key == 'k') {
      eye_move(0,-1,0);
      draw();
    } else if (key == 'l') {
       eye_move(1,0,0);
      draw();
    }else if(key=='u'){
       eye_move(0,0,1);
       draw();
    }else if(key=='o'){
       eye_move(0,0,-1);
       draw();
    }else if(key=='m'){
       eye_rotate(.1,0,0);
       draw();
    }else if(key==','){
       eye_rotate(0,.1,0);
       draw();
    }else if(key=='.'){
       eye_rotate(0,0,.1);
       draw();
    }else if (key > 47 && key < 58) {
      draw();
      k = (key - 48);
    }
  }
  
}
