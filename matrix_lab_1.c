#include <D2d_matrix.h>

  int D2d_print_mat(double m[3][3]) {
    int i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        printf(" %12.4lf ", m[i][j]);
      }
      printf("\n");
    }

    return 1;
  }

int D2d_copy_mat(double a[3][3], double b[3][3]) {
  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      a[i][j] = b[i][j];
    }
  }

  return 1;
}

int D2d_make_identity(double a[3][3]){
    int i, j;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (i == j) a[i][j] = 1.0;
        else a[i][j] = 0.0;
      }
    }

    return 1;
  }

int D2d_mat_mult(double res[3][3], double a[3][3], double b[3][3]) {
  double sum, t[3][3];
  int i, j, k;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      sum = 0;
      for (k = 0; k < 3; k++) {
        sum += (a[i][k] * b[k][j]);
      }
      t[i][j] = sum;
    }
  }

  D2d_copy_mat(res, t);
  return 1;
}

int D2d_mat_mult_points(double * X, double * Y, double m[3][3], double * x, double * y, int n) {
  double temp_x[n];
  int i;
  for (i = 0; i < n; i++) {
    temp_x[i] = (m[0][0] * x[i]) + (m[0][1] * y[i]) + (m[0][2] * 1);
    

  }

  for (i = 0; i < n; i++) {
    Y[i]  = (m[1][0] * x[i]) + (m[1][1] * y[i]) + (m[1][2] * 1);

  }
  
   for (i = 0; i < n; i++) {
    X[i] = temp_x[i];
    }
  return 1;
}

int D2d_scale(double a[3][3], double b[3][3], double xscl, double yscl) {
  double temp[3][3];
  D2d_make_identity(temp);

  temp[0][0] = xscl;
  temp[1][1] = yscl;
  D2d_mat_mult(a, temp, a);

  temp[0][0] = (1 / xscl);
  temp[1][1] = (1 / yscl);
  D2d_mat_mult(b, b, temp);

  return 1;
}

int D2d_rotate(double a[3][3], double b[3][3], double theta) {
  double temp[3][3];
  D2d_make_identity(temp);
  //cos(x) is an even function.
  //sin(x) is an odd function.
  double cosine = cos(theta);
  double sine = sin(theta);
  
  temp[0][0] = cosine;
  temp[0][1] = -1 * sine;
  temp[1][0] = sine;
  temp[1][1] = cosine;
  D2d_mat_mult(a, temp, a);

  //We don't need to recompute the cos(-x) and sin(-x)
  //We can just get it from the even/oddness of the functions.
  temp[0][0] = cosine;
  temp[0][1] = sine;
  temp[1][0] = -sine;
  temp[1][1] = cosine;
  D2d_mat_mult(b, b, temp);

  return 1;
}

int D2d_negate_x(double a[3][3], double b[3][3]) {
  double temp[3][3];
  D2d_make_identity(temp);

  temp[0][0] = -1;
  temp[0][1] = -1;
  D2d_mat_mult(a, temp, a);

  temp[0][0] = -1;
  temp[0][1] = -1;
  D2d_mat_mult(b, b, temp);

  return 1;
}

int D2d_negate_y(double a[3][3], double b[3][3]) {
  double temp[3][3];
  D2d_make_identity(temp);

  temp[1][0] = -1;
  temp[1][1] = -1;
  D2d_mat_mult(a, temp, a);

  temp[1][0] = -1;
  temp[1][1] = -1;
  D2d_mat_mult(b, b, temp);

  return 1;
}

int D2d_translate(double a[3][3], double b[3][3], double dx, double dy){
    double temp[3][3];

    D2d_make_identity(temp);

    temp[0][2] = dx;
    temp[1][2] = dy;
    D2d_mat_mult(a, temp, a);

    temp[0][2] = -dx;
    temp[1][2] = -dy;
    D2d_mat_mult(b, b, temp);

    return 1;
  }
