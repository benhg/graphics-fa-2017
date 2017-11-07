#include "D3d_matrix.h"

int D3d_print_mat (double a[4][4]) {
  int r, c;
  for (r = 0; r < 4; r++ ) {
      for (c = 0; c < 4; c++ ) {
           printf(" %12.4lf ", a[r][c]);//prints 8 digits, 4 after decimal, 12 total
      }
      printf("\n");
  }

  return 1;
} 

int D3d_copy_mat (double a[4][4], double b[4][4]) {
  int r, c;
  for (r = 0; r < 4; r++ ) {
      for (c = 0; c < 4; c++ ) {
           a[r][c] = b[r][c];
      }
  }

  return 1;
} 

int D3d_make_identity (double a[4][4]) {
  int r, c;
  for (r = 0; r < 4; r++ ) {
      for (c = 0; c < 4; c++ ) {
           if (r == c) a[r][c] = 1.0;
               else    a[r][c] = 0.0;
      }
  }

  return 1;
}

int D3d_mat_mult (double res[4][4], double a[4][4], double b[4][4]) {//mayb double check that this works with just changing the four loops
  double t[4][4];
  int r, c, k;
  double tot;
  for (r = 0; r < 4; r++ ) {
    for (c = 0; c < 4; c++ ) {
      tot = 0;
      for (k = 0; k < 4; k++) {
        tot = tot + (a[r][k] * b[k][c]);
      }
      t[r][c] = tot;
    }
  }

  D3d_copy_mat(res, t);
  return 1;
}

int D3d_translate (double a[4][4], double b[4][4], double dx, double dy, double dz) {
  double t[4][4];
  D3d_make_identity(t);

  t[0][3] =  dx;  t[1][3] = dy;  t[2][3] = dz; 
  D3d_mat_mult(a, t, a);

  t[0][3] = -dx;  t[1][3] = -dy;  t[2][3] = -dz; 
  D3d_mat_mult(b, b, t);

  return 1;
}

int D3d_scale (double a[4][4], double b[4][4], double sx, double sy, double sz) {
  double t[4][4];
  D3d_make_identity(t);

  t[0][0] =  sx;  t[1][1] = sy;  t[2][2] = sz;
  D3d_mat_mult(a, t, a);

  t[0][0] = (1/sx);  t[1][1] = (1/sy);  t[2][2] = (1/sz);
  D3d_mat_mult(b, b,t);
  
  return 1;
}

int D3d_rotate_x (double a[4][4], double b[4][4], double radians) {
  double t[4][4];
  D3d_make_identity(t);  

  double cs = cos(radians);
  double sn = sin(radians);
  double negcs = cos(-1 * radians);
  double negsn = sin(-1 * radians);

  t[1][1] =  cs;  t[1][2] = -1 * sn;
  t[2][1] =  sn;  t[2][2] = cs;
  D3d_mat_mult(a, t, a);

  t[1][1] =  negcs;  t[1][2] = -1 * negsn;
  t[2][1] =  negsn;  t[2][2] = negcs;
  D3d_mat_mult(b, b, t);
  
  return 1;
}

int D3d_rotate_y (double a[4][4], double b[4][4], double radians) {
  double t[4][4];
  D3d_make_identity(t);  

  double cs = cos(radians);
  double sn = sin(radians);
  double negcs = cos(-1 * radians);
  double negsn = sin(-1 * radians);

  t[0][0] =  cs;  t[0][2] = sn;
  t[2][0] =  -1 * sn;  t[2][2] = cs;
  D3d_mat_mult(a, t, a);

  t[0][0] =  negcs;  t[0][2] = negsn;
  t[0][2] =  -1 * negsn;  t[2][2] = negcs;
  D3d_mat_mult(b, b, t);
  
  return 1;
}

int D3d_rotate_z (double a[4][4], double b[4][4], double radians) {
  double t[4][4];
  D3d_make_identity(t);  

  double cs = cos(radians);
  double sn = sin(radians);
  double negcs = cos(-1 * radians);
  double negsn = sin(-1 * radians);

  t[0][0] =  cs;  t[0][1] = -1 * sn;
  t[1][0] =  sn;  t[1][1] = cs;
  D3d_mat_mult(a, t, a);

  t[1][1] =  negcs;  t[1][2] = -1 * negsn;
  t[2][1] =  negsn;  t[2][2] = negcs;
  D3d_mat_mult(b, b, t);
  
  return 1;
}

int D3d_cs_rotate_x (double a[4][4], double b[4][4], double cs, double sn) {
  double t[4][4];
  D3d_make_identity(t);  

  t[1][1] =  cs;  t[1][2] = -1 * sn;
  t[2][1] =  sn;  t[2][2] = cs;
  D3d_mat_mult(a, t, a);
/* what do I do about the inverse here??
  t[1][1] =  negcs;  t[1][2] = -1 * negsn;
  t[2][1] =  negsn;  t[2][2] = negcs;
  D3d_mat_mult(b, b, t);
*/  
  return 1;
}

int D3d_cs_rotate_y (double a[4][4], double b[4][4], double cs, double sn) {
  double t[4][4];
  D3d_make_identity(t);  

  t[0][0] =  cs;  t[0][2] = sn;
  t[0][2] =  -1 * sn;  t[2][2] = cs;
  D3d_mat_mult(a, t, a);

/* what do I do about the inverse here??
  t[0][0] =  negcs;  t[0][2] = negsn;
  t[0][2] =  -1 * negsn;  t[2][2] = negcs;
  D3d_mat_mult(b, b, t);
*/ 
  return 1;
}

int D3d_cs_rotate_z (double a[4][4], double b[4][4], double cs, double sn) {
  double t[4][4];
  D3d_make_identity(t);  

  t[0][0] =  cs;  t[0][1] = -1 * sn;
  t[1][0] =  sn;  t[1][1] = cs;
  D3d_mat_mult(a, t, a);

/* what do I do about the inverse here??
  t[1][1] =  negcs;  t[1][2] = -1 * negsn;
  t[2][1] =  negsn;  t[2][2] = negcs;
  D3d_mat_mult(b, b, t);
*/  
  return 1;
}

int D3d_negate_x (double a[4][4], double b[4][4]) {
  double t[4][4] ;
  D3d_make_identity(t) ; 
  
  t[0][0] = -1;
  D3d_mat_mult(a,  t,a);

  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}

int D3d_negate_y (double a[4][4], double b[4][4]) {
  double t[4][4] ;
  D3d_make_identity(t) ; 
  
  t[1][1] =  -1;
  D3d_mat_mult(a,  t,a);

  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}

int D3d_negate_z (double a[4][4], double b[4][4]) {
  double t[4][4] ;
  D3d_make_identity(t) ; 
  
  t[2][2] =  -1;
  D3d_mat_mult(a,  t,a);

  D3d_mat_mult(b,  b,t) ;

  return 1 ;
}

int D3d_mat_mult_points (double *X, double *Y, double *Z, double m[4][4],
                         double *x, double *y, double *z, int numpoints) {
  double tx[numpoints]; double ty[numpoints]; double tz[numpoints];
  int c;
  double tot;
  for (c = 0; c < numpoints; c++ ) {
    tot = 0;
    tot = (m[0][0] * x[c]);
    tot = tot + (m[0][1] * y[c]);
    tot = tot + (m[0][2] * z[c]);
    tot = tot + (m[0][3] * 1);
    tx[c] = tot;
  }

  for (c = 0 ; c < numpoints ; c++ ) {
    tot = 0;
    tot = (m[1][0] * x[c]);
    tot = tot + (m[1][1] * y[c]);
    tot = tot + (m[1][2] * z[c]);
    tot = tot + (m[1][3] * 1);
    ty[c] = tot;
  }
  
  for (c = 0 ; c < numpoints ; c++ ) {
    tot = 0;
    tot = (m[2][0] * x[c]);
    tot = tot + (m[2][1] * y[c]);
    tot = tot + (m[2][2] * z[c]);
    tot = tot + (m[2][3] * 1);
    tz[c] = tot;
  }

  for (c = 0; c < numpoints; c++) {
    X[c] = tx[c];
    Y[c] = ty[c];
    Z[c] = tz[c];
  }

  return 1;
}

int D3d_mat_mult_pt (double P[3], double m[4][4], double Q[3]) {
  double tp[3];
  int c;
  double tot;
  
  for(c= 0; c < 3; c++) {
    tot = 0;
    tot = (m[c][0] * Q[0]);
    tot = tot + (m[c][1] * Q[1]);
    tot = tot + (m[c][2] * Q[2]);
    tot = tot + (m[c][3] * 1);
    tp[c] = tot;
  }  
  
  for (c = 0; c < 3; c++) {
    P[c] = tp[c];
  }  
}

int D3d_x_product (double res[3], double a[3], double b[3]) {
  res[0] = ((a[1] * b[2]) - (b[1] * a[2]));
  res[1] = (-1 * ((a[0] * b[2]) - (b[0] * a[2])));
  res[2] = ((a[0] * b[1]) - (b[0] * a[1]));
}




























