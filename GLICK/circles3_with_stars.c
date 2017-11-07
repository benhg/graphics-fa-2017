#include <FPT.h>

void star (double xc, double yc, double r)
{
  double x[5],y[5],t ;
  int k ;

  for (k = 0 ; k < 5 ; k++) {
    t = k * 2*2*M_PI/5 + M_PI/2 ; 
    x[k] = xc + r*cos(t) ;
    y[k] = yc + r*sin(t) ;
  }
  
  G_polygon(x,y,5) ;
}



void fill_star_flaw (double xc, double yc, double r)
{
  double x[5],y[5],t ;
  int k ;

  for (k = 0 ; k < 5 ; k++) {
    t = k * 2*2*M_PI/5 + M_PI/2 ; 
    x[k] = xc + r*cos(t) ;
    y[k] = yc + r*sin(t) ;
  }
  
  G_fill_polygon(x,y,5) ; // fun but not correct
}



void fill_star (double xc, double yc, double r)
{
  double x[5],y[5],t ;
  int k ;

  for (k = 0 ; k < 5 ; k++) {
    t = k * 2*2*M_PI/5 + M_PI/2 ; 
    x[k] = xc + r*cos(t) ;
    y[k] = yc + r*sin(t) ;
  }
  
  int j = 1 ;
  for (k = 0 ; k < 5 ; k++) {
    G_fill_triangle(xc,yc, x[k],y[k], x[j],y[j]) ;
    j = (j + 1) % 5 ;
  }

}



int main()
{
  double a[2],b[2],dx,dy,r ;
  G_init_graphics(600,600) ;

  G_rgb(0,0,0) ;
  G_clear() ;

  G_rgb(1,1,0) ;
  G_fill_rectangle(0,0,20,20) ;

  G_rgb(1,0,1) ;
  while (0 == 0) {
    G_wait_click(a) ;   G_fill_circle(a[0],a[1],2) ;
    if ((a[0] < 20) && (a[1] < 20)) { break ; }
    G_wait_click(b) ;   // G_fill_circle(b[0],b[1],2) ;
    dx = b[0] - a[0] ; dy = b[1] - a[1] ; r = sqrt(dx*dx + dy*dy) ;
    G_circle(a[0],a[1],r) ;
    star(a[0],a[1],r) ;
    //   fill_star(a[0],a[1],r) ;
  }

  //  G_wait_key() ;
}
