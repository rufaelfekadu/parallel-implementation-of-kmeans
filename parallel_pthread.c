#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#define MATCH(s) (!strcmp(argv[ac], (s)))

int max_range = 800;
int max_iterations = 100;
pthread_mutex_t mutexcluster;


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

struct thread_data
{
    clusters * cl;
    points * pt;
    int tid;
    int *num_thread;
    
};


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

    // double min_distance;
    // int min_index;

    // for(int i = 0; i < points_size; i++){

    //     min_distance = DBL_MAX;
    //     min_index = 0;

    //     for(int j = 0; j < clusters_size; j++){


    //         double distance = euclidean_dist(&po->p[i], &clust->c[j]);

    //         if(distance < min_distance){

    //             min_distance = distance;
    //             min_index = j;
    //         }
    //     }
    //     po->p[i].cluster_id = min_index;
    //     // points[i].set_cluster_id(min_index);
    //     #pragma omp atomic
    //     clust->c[min_index].new_x_coord += po->p[i].x_coord;
    //     #pragma omp atomic
    //     clust->c[min_index].new_y_coord += po->p[i].y_coord;
    //     #pragma omp atomic
    //     clust->c[min_index].size++;
    //     // clusters[min_index].add_point(points[i]);
    // }


}

void * compute_thread(void*  thread_data){

    struct thread_data *data = (struct thread_data *) thread_data; 
    unsigned long points_size = data->pt->size;
    unsigned long clusters_size = data->cl->size;
    int* num = data->num_thread;
    int t = data->tid;



    double min_distance;
    int min_index;
    int start = t * (points_size/ *num);
    int end = (t+1)*(points_size/ *num);
    if (end== *num-1)
    {
        end = points_size;
    }

    for(int i = start; i < end; i++){

        min_distance = DBL_MAX;
        min_index = 0;
        
        for(int j = 0; j < clusters_size; j++){


            double distance = euclidean_dist(&data->pt->p[i], &data->cl->c[j]);

            if(distance < min_distance){

                min_distance = distance;
                min_index = j;
            }
        }
        data->pt->p[i].cluster_id = min_index;
        // points[i].set_cluster_id(min_index);
        // #pragma omp atomic
        pthread_mutex_lock(&mutexcluster);
        data->cl->c[min_index].new_x_coord += data->pt->p[i].x_coord;
        // #pragma omp atomic
        data->cl->c[min_index].new_y_coord += data->pt->p[i].y_coord;
        // #pragma omp atomic
        data->cl->c[min_index].size++;
        pthread_mutex_unlock(&mutexcluster);
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

    for(int i = 0; i < poin->size; i++){

        fprintf(fout, "%d %d %d\n",poin->p[i].x_coord,poin->p[i].y_coord,poin->p[i].cluster_id);

    }

    fclose(fout);
    system("gnuplot -p -e \"plot 'data.txt' using 1:2:3 with points palette notitle\"");
    remove("data.txt");

}

void result(double dur,int size,int num_cluster){
    char dir[105]="./output/result.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",num_cluster, dur);

    fclose(fout);
}

// void pthread_vs_cluster(double dur,int num_cluster){
//     char dir[105]="./output/to_plot/pthread_vs_cluster.txt";
// 	FILE *fout = fopen(dir, "a");
// 	fprintf(fout, "%d %f\n",num_cluster, dur);

//     fclose(fout);
// }

void pthread_vs_point(double dur,int size){
    char dir[105]="./output/to_plot/pthread_vs_point.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",size, dur);

    fclose(fout);
}

void pthread_vs_thread(double dur,int numthread){
    char dir[105]="./output/to_plot/pthread_vs_thread.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",numthread, dur);

    fclose(fout);
}

void plot_result(){

    system("gnuplot -p -e \"plot 'output/result.txt' using 1:2:3 with points palette notitle\"");

}
int main(int argc,char **argv){

    pthread_mutex_init(&mutexcluster, NULL);
    srand(31359);
    int ac,numthreads = 1;
    int disable_display = 0;
    int seedVal = 100;
    int nx = 1;
    int num_cluster = 4;

    for(ac=1;ac<argc;ac++)
    {
        if(MATCH("-n")) {nx = atoi(argv[++ac]);}
        else if(MATCH("-i")) {max_iterations = atoi(argv[++ac]);}
        else if(MATCH("-t"))  {numthreads = atof(argv[++ac]);}
        else if(MATCH("-c"))  {num_cluster = atof(argv[++ac]);}
        // else if(MATCH("-s"))  {seedVal = atof(argv[++ac]);}
        else if(MATCH("-d"))  {disable_display = 1;}
        else {
            printf("Usage: %s [-n <dataset>] [-i <iterations>] [-t numthreads] [-c num clusters] [-d disable display]\n",argv[0]);
            return(-1);
        }
    }

    char *dataset_filename;
    switch (nx)
    {
    case 1:
        dataset_filename = "dataset-10000.txt";
        break;
    case 2:
        dataset_filename = "dataset-50000.txt";
        break;
    case 3:
        dataset_filename = "dataset-100000.txt";
        break;
    case 4:
        dataset_filename = "dataset-200000.txt";
        break;
    case 5:
        dataset_filename = "dataset-400000.txt";
        break;
    case 6:
        dataset_filename = "dataset-500000.txt";
        break;
    case 7:
        dataset_filename = "dataset-600000.txt";
        break;
    case 8:
        dataset_filename = "dataset-800000.txt";
        break;
    case 9:
        dataset_filename = "dataset-1000000.txt";
        break;
    default:
        dataset_filename = "dataset-10000.txt";
        break;
    }



    double time_point1 = omp_get_wtime();

    printf("Starting initialization..\n");

    printf("Creating points..\n");
    points mypoints = init_points(dataset_filename,100);
    printf("Points initialized \n");

    
    printf("Creating clusters..\n");
    clusters mycluster = init_clusters(num_cluster);
    printf("Clusters initialized \n");

    double time_point2 = omp_get_wtime();
    double duration = time_point2 - time_point1;

    printf("Points and clusters generated in: %f seconds\n", duration);
    printf("number of points: %d\n",mypoints.size);
    printf("number of cluster: %d\n",mycluster.size);

    bool conv = true;
    int iterations = 0;

    printf("Starting iterate..\n");
    double time_point3 = omp_get_wtime();
    while(conv && iterations < max_iterations){

        iterations ++;

        // compute_distance(&mypoints, &mycluster);
        pthread_t threads[numthreads];
        struct thread_data thread_data_array[numthreads];

        for(int t=0; t<numthreads; t++){
            thread_data_array[t].cl = &mycluster;
            thread_data_array[t].pt = &mypoints;
            thread_data_array[t].tid = t;
            thread_data_array[t].num_thread = &numthreads;
            int rc = pthread_create(&threads[t], NULL, compute_thread, (void *) &thread_data_array[t]);
        }

        for(int i = 0; i < numthreads; i++)
        {
            if(pthread_join(threads[i], NULL))
            {
                    printf("Could not join thread %d\n", i);
                    return -1;
            }
        }

        // pthread_exit(NULL);
        conv = update_clusters(&mycluster);

        // printf("Iteration %d done \n", iterations);

    }

    double time_point4 = omp_get_wtime();
    duration = time_point4 - time_point3;

    printf("Number of iterations: %d, total time: %f seconds, time per iteration: %f seconds\n",
        iterations, duration, duration/iterations);

    // result(duration/iterations,mypoints.size,num_cluster);
    pthread_vs_point(duration/iteration,mypoints.size);
    pthread_vs_cluster(duration/iteration,num_cluster);
    pthread_vs_thread(duration/iteration,numthreads);

    if(!disable_display){

        printf("Drawing the chart...\n");
        draw_chart_gnu(&mypoints);

    }
    

    free(mypoints.p);
    free(mycluster.c);

    return 0;
}