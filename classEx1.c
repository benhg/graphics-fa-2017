#include <FPT.h>
#include <math.h>

int main(){
  double p[2], r[2]:
  G_init_graphics(500, 500);
  G_wait_click(p);
  G_wait_click(r);
  double rad = sqrt((r[0]-p[0])^2 + (r[1]-p[1])^2)
  G_circle(p[0], p[1], rad)
  G_wait_key()
  G_close()


}
