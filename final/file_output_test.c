#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main()
{
  int n ;
  double x ;
  FILE *fp ;

  fp = fopen("test.xyz", "w") ;  // the "w" indicates we wish to write the file
  if (fp == NULL) {
    printf("can't open the file\n") ;
    exit(1) ;
  }
  
  n = 12 ;
  x = 4.9 ;

  fprintf(fp,  "%d  %lf\n",n,x) ;

}
