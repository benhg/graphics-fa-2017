#include <FPT.h>
#include <D2d_matrix.h>

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


int main()
{


  G_init_graphics(600,600) ;

  G_rgb(0,0,0) ;  G_clear() ;
  G_rgb(0,0,1) ;  G_fill_polygon(jx,jy,5) ;
  G_rgb(1,1,0) ;  G_draw_string("any key to continue", 250,10) ;
  G_rgb(1,0,1) ;  G_line(0,500,500,0) ;
  G_wait_key() ;


}
