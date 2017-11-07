#include <FPT.h>
#include "D2d_matrix.h"
#define HEIGHT 600
#define WIDTH 600
// Run the program, jetm.exe, by typing  ./jetm.exe

// By continually holding down on a key, you will see a movie.
// Your task is to duplicate this program.

// Notice that the jet keeps moving until the user types the 'q' key,
// and then the program terminates.  Also notice that the program
// prints the number of each frame that is displayed.

// With each key-press, the jet moves 2 units (each pixel is 1 unit by 1 unit).

// You must move the jet by only using matrix operations-NOT procedurally.

// For full credit, compute any 3x3 matrices you need BEFORE the loop
// that displays the movie.  An efficient solution does not need
// to perpetually rebuild them with each new scene.


//jet
double jx[5] = {0, 40, 35, 10,  0} ;
double jy[5] = {0,  0,  8,  8, 15} ;
double kevin[3][3];

void stretch(double sx, double sy) {
  double m[3][3], minv[3][3];
  D2d_make_identity(m);
  D2d_make_identity(minv);
  D2d_scale(m, minv,  sx, sy);
  D2d_mat_mult_points(jx, jy,  m, jx, jy,5);
}

void translate(double dx, double dy) {
  double m[3][3], minv[3][3];
  D2d_make_identity(m);
  D2d_make_identity(minv);
  D2d_translate(m, minv, dx, dy);
  D2d_mat_mult_points(jx, jy,  m, jx, jy, 5);
}

void rotate(double radians) {
  double m[3][3], minv[3][3];
  D2d_make_identity(m); D2d_make_identity(minv);
  D2d_rotate(m, minv,  radians);
  D2d_mat_mult_points(jx, jy,  m, jx, jy,5);
}

int D2d_flip(double a[3][3]) {
  double temp[3][3];
  D2d_make_identity(temp);
  temp[0][0] = 1;
  temp[0][1] = 0;
  temp[1][0] = 0;
  temp[1][1] = -1;
  D2d_mat_mult(a, temp, a);
    return 1;
}

void flip(){
   D2d_mat_mult_points(jx, jy,  kevin, jx, jy, 5);
}

int main()
{
  

  
  G_init_graphics(600,600) ;

  G_rgb(0,0,0) ;  G_clear() ;
  G_rgb(0,0,1) ;
  G_rgb(1,1,0) ;  G_draw_string("any key to continue", 250,10) ;
  G_rgb(1,0,1) ;  G_line(0,500,500,0) ;
  G_wait_key() ;
  rotate(-.8);
  translate(0,500);
   D2d_make_identity(kevin);
   D2d_flip(kevin);
  int i, key;
  i = 0;
  while(1){
    key = G_wait_key();
    if(key =='q'){
      break ;
    }
    G_rgb(0,0,0);
    G_clear();
    G_rgb(1,0,1);
    G_line(0,500,500,0) ;
    G_rgb(0,0,1);
    G_fill_polygon(jx,jy,5) ;
    translate(2,-2);
    if(i%50==0 && i>0){
      rotate(.8);
      double cheight;
      cheight = jy[0];
      translate(0,-cheight);
      flip();
      translate(0, cheight);
      rotate(-.8);
    }
    printf("%d\n", i);
    i++;
  }

}
