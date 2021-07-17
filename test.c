#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <math.h>
#include <time.h>


int max_range = 800;
int max_iterations = 100;

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

typedef enum {false, true} bool;

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
    for(int i=0;i<po.size;i++){
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

        cl.c[i] = init_cluster((rand()-RAND_MAX/2) % (int) max_range, (rand()-RAND_MAX/2) % (int) max_range);
        cl.items++;

    }

   return cl; 
}

double euclidean_dist(point* point, cluster* clust){

    double distance = sqrt(pow(point->x_coord - clust->x_coord,2) +
                           pow(point->y_coord - clust->y_coord, 2));

    return distance;
}

void compute_distance(points* po, clusters* clust){

    // printf("%f\n",clust->c[4].x_coord);
    unsigned long points_size = po->size;
    unsigned long clusters_size = clust->size;

    double min_distance;
    int min_index;

    for(int i = 0; i < points_size; i++){

        min_distance = DBL_MAX;
        min_index = 0;

        for(int j = 0; j < clusters_size; j++){


            double distance = euclidean_dist(&po->p[i], &clust->c[j]);

            if(distance < min_distance){

                min_distance = distance;
                min_index = j;
            }
        }
        po->p[i].cluster_id = min_index;
        // points[i].set_cluster_id(min_index);
        #pragma omp atomic
        clust->c[min_index].new_x_coord += po->p[i].x_coord;
        #pragma omp atomic
        clust->c[min_index].new_y_coord += po->p[i].y_coord;
        #pragma omp atomic
        clust->c[min_index].size++;
        // clusters[min_index].add_point(points[i]);
    }
}

void free_point(cluster* clu){
    clu->size = 0;
    clu->new_x_coord = 0;
    clu->new_y_coord = 0;
}

bool update_coords(cluster * clus){


    if(clus->x_coord == clus->new_x_coord/clus->size && clus->y_coord == clus->new_y_coord/clus->size){
        return false;
    }

    clus->x_coord = clus->new_x_coord/clus->size;
    clus->y_coord = clus->new_y_coord/clus->size;

    return true;

}

bool update_clusters(clusters *clus){

    bool conv = false;

    for(int i = 0; i < clus->size; i++){

        conv = update_coords(&clus->c[i]);
        free_point(&clus->c[i]);
    }

    return conv;
}

//Draw point plot with gnuplot
void draw_chart_gnu(points* poin){

    // char dir[105] = "./output/";
    // strcat(dir,cluster_filename);
	FILE *fout = fopen("data.txt", "w");

    // ofstream outfile("data.txt");

    for(int i = 0; i < poin->size; i++){

        fprintf(fout, "%d %d %d\n",poin->p[i].x_coord,poin->p[i].y_coord,poin->p[i].cluster_id);

        // outfile << point.get_x_coord() << " " << point.get_y_coord() << " " << point.get_cluster_id() << std::endl;

    }

    fclose(fout);
    // outfile.close();
    system("gnuplot -p -e \"plot 'data.txt' using 1:2:3 with points palette notitle\"");
    // remove("data.txt");

}

void result(double dur,int size){
    char dir[105]="./output/result.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f 1\n", size, dur);

    fclose(fout);
}

void plot_result(){

    system("gnuplot -p -e \"plot 'output/result.txt' using 1:2:3 with points palette notitle\"");

}
int main(){

    // points mypoints = init_points("dataset-10000.txt",100);
    // clusters mycluster = init_clusters(4);
    // int i=0;
    // while(i<mycluster.size){
    //     printf("%f %f %d\n",mycluster.c[i].x_coord,mycluster.c[i].y_coord,i);
    //     i++;
    // }

    // printf("Number of points %d\n", num_point);
    // printf("Number of clusters %d\n", num_cluster);

    srand(time(NULL));

    double time_point1 = omp_get_wtime();

    printf("Starting initialization..\n");

    printf("Creating points..\n");
    points mypoints = init_points("dataset-1000000.txt",100);
    printf("Points initialized \n");

    
    printf("Creating clusters..\n");
    clusters mycluster = init_clusters(4);
    printf("Clusters initialized \n");

    double time_point2 = omp_get_wtime();
    double duration = time_point2 - time_point1;

    printf("Points and clusters generated in: %f seconds\n", duration);
    printf("number of points: %d\n",mypoints.size);
    printf("number of cluster: %f\n",mycluster.c[2].x_coord);

    bool conv = true;
    int iterations = 0;

    printf("Starting iterate..\n");

     while(conv && iterations < max_iterations){

        iterations ++;

        compute_distance(&mypoints, &mycluster);

        conv = update_clusters(&mycluster);

        // printf("Iteration %d done \n", iterations);

    }

    double time_point3 = omp_get_wtime();
    duration = time_point3 - time_point2;

    printf("Number of iterations: %d, total time: %f seconds, time per iteration: %f seconds\n",
           iterations, duration, duration/iterations);

    result(duration,mypoints.size);
    printf("Drawing the chart...\n");
    draw_chart_gnu(&mypoints);
    plot_result();

    // free(&mypoints);
    // free(&mycluster);
    return 0;
}