#include <stdio.h>
#include <stdlib.h>
// #include <mpi.h>
#include <assert.h>


float* init_points(const int num_elements, int site_per_proc,int nprocs) {
  float *X,*Y;
  X = (float *)malloc(sizeof(float) * num_elements);
  Y = (float *)malloc(sizeof(float) * num_elements);
  assert(X != NULL);
  assert(Y != NULL);
  float *all_points = (float *)malloc(sizeof(float) * num_elements * 2);
    int * num;
  FILE *fin = fopen("./datasets/dataset-10000.txt", "r");
	fscanf(fin, "%d", num);
    for(int i=0;i<num_elements;i++){
      
        fscanf(fin, "%f %f",&X[i],&Y[i]);
    }
  fclose(fin);
  int m=0;
  int x=0;
  int y=0;
  for (int j=0;j<nprocs;j++ ){

    for(int k=0;k<site_per_proc;k++){
      all_points[m]=X[x];
      x++;
      m++;
    }
    for(int k=0;k<site_per_proc;k++){
      all_points[m]=Y[y];
      y++;
      m++;
    }

  }
  free(X);
  free(Y);
  return all_points;
}
int main(int argc, char const *argv[])
{
    float * a = NULL;
    // FILE *fout = fopen("output/data2.txt", "w");
    a = init_points(10000,100,100);
    int i;
    for (i=0;i<100;i++){
        printf("%f ",a[i]);
        // fprintf(fout, "%4f\n",a[i]);
    }
    printf("%d",i);
    // fclose(fout);

    return 0;
}
