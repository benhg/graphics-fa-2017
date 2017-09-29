#include <fpt.h>
int n;
double x[9001], y[9001];
int psize[5000];
int con[5000][50];
double red[5000];
double green[5000];
double blue[5000];

int main(){
  FILE *q;
  int i;
  q = fopen("argv");
  if(q==NULL){
    printf("Can't open file\n");
    exit(1);
  }
  fscanf(q, "%d", &numpoints);
  for(i=0;i<numpoints;i++){
    fscanf(q,"%lf %lf",&x[i],&y[i]);
  }
}
