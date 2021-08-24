#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define MATCH(s) (!strcmp(argv[ac], (s)))

int max_range = 800;
int max_iterations = 100;
int numberOfElements;
int num_process;
int numOfClusters = 4;

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
    numberOfElements = po.size;
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

    for(int i = 0; i < poin->size; i++){

        fprintf(fout, "%d %d %d\n",poin->p[i].x_coord,poin->p[i].y_coord,poin->p[i].cluster_id);

    }

    fclose(fout);
    system("gnuplot -p -e \"plot 'data.txt' using 1:2:3 with points palette notitle\"");
    // remove("data.txt");

}

void result(double dur,int size,int num_cluster){
    char dir[105]="./output/result.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",num_cluster, dur);

    fclose(fout);
}

void plot_result(){

    system("gnuplot -p -e \"plot 'output/result.txt' using 1:2:3 with points palette notitle\"");

}
int main(int argc,char **argv){

    // initialize the MPI environment
	MPI_Init(NULL, NULL);

	// get number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// get rank
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    points mypoints;
    clusters myclusters;

    // send buffers
	double *k_means_x = NULL;		// k means corresponding x values
	double *k_means_y = NULL;		// k means corresponding y values
	int *k_assignment = NULL;		// each data point is assigned to a cluster

	// receive buffer
	int *recv_x = NULL;
	int *recv_y = NULL;
	int *recv_assign = NULL;

    if(world_rank==0){
        srand(50);
        int ac,numthreads = 1;
        int disable_display = 0;
        int seedVal = 100;
        int nx = 1;

        for(ac=1;ac<argc;ac++)
        {
            if(MATCH("-n")) {nx = atoi(argv[++ac]);}
            else if(MATCH("-i")) {max_iterations = atoi(argv[++ac]);}
            else if(MATCH("-p")) {num_process = atoi(argv[++ac]);}
            else if(MATCH("-t"))  {numthreads = atof(argv[++ac]);}
            else if(MATCH("-c"))  {numOfClusters = atof(argv[++ac]);}
            else if(MATCH("-s"))  {seedVal = atof(argv[++ac]);}
            else if(MATCH("-d"))  {disable_display = 1;}
            else {
                printf("Usage: %s [-n < meshpoints>] [-i <iterations>] [-t numthreads] [-s seed] [-p prob] [-d]\n",argv[0]);
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

        MPI_Bcast(&numOfClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);

        printf("Starting initialization..\n");

        printf("Creating points..\n");
        mypoints = init_points(dataset_filename,100);
        printf("Points initialized \n");


        printf("There are a total number of %d elements in the file.\n", numberOfElements);

		// broadcast the number of elements to all nodes
		MPI_Bcast(&numberOfElements, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        printf("Creating clusters..\n");
        myclusters = init_clusters(4);
        printf("Clusters initialized \n");

        // printf("Points and clusters generated in: %f seconds\n", duration);
        // printf("number of points: %d\n",mypoints.size);
        // printf("number of cluster: %d\n",mycluster.size);

        printf("Starting iterate..\n");
        double time_point3 = omp_get_wtime();

        // allocate memory for receive buffers
		recv_x = (int *)malloc(sizeof(int) * ((numberOfElements/num_process) + 1));
		recv_y = (int *)malloc(sizeof(int) * ((numberOfElements/num_process) + 1));
		recv_assign = (int *)malloc(sizeof(int) * ((numberOfElements/num_process) + 1));

        if(recv_x == NULL || recv_y == NULL || recv_assign == NULL)
		{
			perror("malloc");
			exit(-1);
		}

    }else{	// I am a worker node
        int ac;
        for(ac=1;ac<argc;ac++)
        {
            if(MATCH("-i")) {max_iterations = atoi(argv[++ac]);}
            else if(MATCH("-p")) {num_process = atoi(argv[++ac]);}
            else {
                printf("Usage: %s [-n < meshpoints>] [-i <iterations>] [-t numthreads] [-s seed] [-p prob] [-d]\n",argv[0]);
                return(-1);
            }
        }

		// num_process = atoi(argv[1]);

		// receive broadcast of number of clusters
		MPI_Bcast(&numOfClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// receive broadcast of number of elements
		MPI_Bcast(&numberOfElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for arrays
		// k_means_x = (double *)malloc(sizeof(double) * numOfClusters);
		// k_means_y = (double *)malloc(sizeof(double) * numOfClusters);

		// if(k_means_x == NULL || k_means_y == NULL)
		// {
		// 	perror("malloc");
		// 	exit(-1);
		// }

		// allocate memory for receive buffers
		recv_x = (int *)malloc(sizeof(int) * ((numberOfElements/num_process) + 1));
		recv_y = (int *)malloc(sizeof(int) * ((numberOfElements/num_process) + 1));
		recv_assign = (int *)malloc(sizeof(int) * ((numberOfElements/num_process) + 1));

		if(recv_x == NULL || recv_y == NULL || recv_assign == NULL)
		{
			perror("malloc");
			exit(-1);
		}
	}

    int *x;
    int *y;
    for(int i=0;i < mypoints.size;i++){
        x[i] = &mypoints.p[i].x_coord;
        y[i] = &mypoints.p[i].y_coord; 
    }
    MPI_Scatter(x, (numberOfElements/num_process) + 1, MPI_INT,
		recv_x, (numberOfElements/num_process) + 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Scatter(x, (numberOfElements/num_process) + 1, MPI_INT,
		recv_x, (numberOfElements/num_process) + 1, MPI_INT, 0, MPI_COMM_WORLD);

	// MPI_Scatter(data_y_points, (numOfElements/num_of_processes) + 1, MPI_DOUBLE,
	// 	recv_y, (numOfElements/num_of_processes) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    bool conv = true;
    int iterations = 0;
    double time_point3 = omp_get_wtime();

    while(conv && iterations < max_iterations){

        MPI_Bcast(k_means_x, numOfClusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(k_means_y, numOfClusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        iterations ++;

        compute_distance(&mypoints, &myclusters);

        conv = update_clusters(&myclusters);

        // printf("Iteration %d done \n", iterations);

    }

    double time_point4 = omp_get_wtime();
    duration = time_point4 - time_point3;

    printf("Number of iterations: %d, total time: %f seconds, time per iteration: %f seconds\n",
        iterations, duration, duration/iterations);

    result(duration/iterations,mypoints.size,num_cluster);

    if(!disable_display){

        printf("Drawing the chart...\n");
        draw_chart_gnu(&mypoints);

    }
    

    free(mypoints.p);
    free(mycluster.c);

    return 0;
}