#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int max_range = 1000000;

struct Cluster{

    double x_coord;
    double y_coord;
    //Accumulate the Point coords here
    double new_x_coord;
    double new_y_coord;
    //Number of points inside this Cluster
    int size;

};
typedef struct Cluster cluster;

struct Point{

  int x_coord;
  int y_coord; 
  int cluster_id; 

};
typedef struct Point point;

struct Points
{
    point* p;
    int size;
    int items;
};
typedef struct Points points;

struct Clusters{
    cluster* c;
    int size;
    int items;
};
typedef struct Clusters clusters;

point init_point(int x_cord,int y_cord){
    point p;
    p.cluster_id=0;
    p.x_coord = x_cord;
    p.y_coord = y_cord;
return p;
}

cluster init_cluster(int x_cord,int y_cord){

    cluster  c;
    c.new_x_coord = 0;
    c.new_y_coord = 0;
    c.size = 0;
    c.x_coord = x_cord;
    c.y_coord = y_cord;

    return c;
}
points init_points(const char *dataset_filename,int num_points){
    points po;
    char dir[105]="./datasets/";
    strcat(dir,dataset_filename); 
	FILE *fin = fopen(dir, "r");
	fscanf(fin, "%d", &po.size);

    po.p = (point*)malloc(sizeof(point)*po.size);
    po.items =0;
    for(int i=0;i<num_points;i++){
        int x,y;
        fscanf(fin, "%d %d",&x,&y);
        po.p[i] = init_point(x,y);
        po.items++;
    }

    return po;
}

clusters init_clusters(int num_cluster){

    clusters cl;
    cl.size=num_cluster;
    cl.c = (cluster*)malloc(sizeof(cluster)*cl.size);
    cl.items=0;
    for(int i = 0; i < cl.size; i++){

        cl.c[i] = init_cluster(rand() % (int) max_range, rand() % (int) max_range);
        cl.items++;

    }

   return cl; 
}

int main(){

    points mypoints = init_points("dataset-10000.txt",100);
    clusters mycluster = init_clusters(4);
    int i=0;
    while(i<mycluster.size){
        printf("%f %f %d\n",mycluster.c[i].x_coord,mycluster.c[i].y_coord,i);
        i++;
    }

    return 0;
}