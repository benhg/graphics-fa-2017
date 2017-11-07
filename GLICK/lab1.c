
#include <FPT.h>
#include <math.h>
#include <stdlib.h>
#define GRAPHICS_WIDTH 1000
#define GRAPHICS_HEIGHT 1000
#define TRUE 1
#define FALSE !TRUE

void polygon(double *x, double *y, int n){
	int i, j;
	for(i = 0; i < n ; i++){
	j = (i+1) % n;
	G_line(x[i],y[i], x[j], y[j]);	
	}
}
double intersect(double *a, double *b, double y){
	if(a[1] < y && b[1] > y || a[1] > y && b[1]< y){
	double m = (b[1]-a[1])/(b[0]-a[0]);
	double off = b[1] - (b[0]*m);
	return (y-off) / m;
	}
}

double* sort_double(double *grade, int n){
  int i, j, swapped, k;
  double temp;
  for (i = 0; i < n; ++i)
    {
      for (j = i + 1; j < n; ++j)
	{
	  if (grade[i] < grade[j])
	    {
	      temp = grade[i];
	      grade[i] = grade[j];
	      grade[j] = temp;
	    }//end if
	}//end inner for
    }//end outer for

  return grade;
}//end sortGrade

void fill_polygon(double *x, double *y, int n){
  int j, k,i;
	for (i=0; i < GRAPHICS_HEIGHT; i++){
	  double points[n];
	  double point1[2], point2[2];
	  for(j=0; j < n; j++){
	    point1[0] = x[j];
	    point1[1] = y[j];
	    int val0;
	    val0 = (j+1)%n;
	    point2[0] = x[val0];
	    point2[1] = y[val0];
	    double l;
	    l=i+.1;
	    points[j] = intersect(point1, point2, l);
	    
	  }
	  double* sorted;
	  sorted = sort_double(points, n);
	  //rmdup(sorted, n);
	  int l;
	  for(l=0;l<n;l+=2){
	    G_line(sorted[l], i, sorted[l+1], i);
	  }
	  
	}
}

void draw(double x_c, double y_c, double rad){
  G_circle(x_c, y_c, rad);
  double x[10],y[10];
  double angle = 1.26;
  double theta = 3.1415 / 2;
  int i = 0;
  for(theta = 3.1415 / 2.0; theta < (3.1415 / 2 ) + (8 * angle); theta += 2 * angle){
    printf("%f : %f\n", theta, rad);
    x[i] = (cos(theta) * rad) + x_c;
    y[i] = (sin(theta) * rad) + y_c;
    printf("%f,%f\n", x[i], y[i]);
    i++;
  }
  // circle(x[0], y[0], 10);
  fill_polygon(x, y, 5);
}

int main(){
  G_init_graphics(GRAPHICS_WIDTH, GRAPHICS_HEIGHT);
  G_rgb(0, 0, 0);
  int num_sides;
  num_sides=0;
  while(TRUE){
    double p[2], r[2];
    double p_x[100];
    double p_y[100];
    G_fill_rectangle(0,0,50,50);
    G_fill_rectangle(GRAPHICS_WIDTH-50, 0, GRAPHICS_WIDTH, 50);
    G_fill_rectangle(GRAPHICS_WIDTH-50, GRAPHICS_HEIGHT-50, GRAPHICS_WIDTH, GRAPHICS_HEIGHT);
    G_fill_rectangle(0, GRAPHICS_HEIGHT-50, 50, GRAPHICS_HEIGHT);
    G_wait_click(p);
    if(p[0] < 50 && p[1] < 50){
      break;
    }
    if(p[0] > GRAPHICS_WIDTH-50 && p[1] >  GRAPHICS_HEIGHT-50){
      G_rgb(0,0,255);
      fill_polygon(p_x, p_y, num_sides);
      G_rgb(0,0,0);
      num_sides = 0;
    }
    if((p[0] < 50) && ( p[1] >  GRAPHICS_HEIGHT-50)){
      G_rgb(255,0,0);
      G_polygon(p_x, p_y, num_sides);
      G_rgb(0,0,0);
    }
    
    if(!(p[0] > GRAPHICS_WIDTH-50) && !( p[1] >  GRAPHICS_HEIGHT-50) && !(p[0] > GRAPHICS_WIDTH-50) && !( p[1] < 50)){
      p_x[num_sides] = p[0];
      p_y[num_sides] = p[1];
      num_sides+=1;
    }
     if(!(p[0] < GRAPHICS_WIDTH-50) && !( p[1] > 50)){
       G_rgb(255,255,255);
       G_clear();
       G_rgb(0,0,0);
    }
	
    G_point(p[0], p[1]);
  }
  G_close();
	
}
