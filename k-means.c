/*  This is an implementation of the k-means clustering algorithm (aka Lloyd's algorithm) using MPI (message passing interface). */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>
#include<math.h>
#include<errno.h>
#include<mpi.h>
#include<time.h>

#define MAX_ITERATIONS 1000

int numOfClusters = 0;
int numOfElements = 0;
int num_of_processes = 0;
int max_range =800;

/* This function goes through that data points and assigns them to a cluster */
void assign2Cluster(double k_x[], double k_y[], double recv_x[], double recv_y[], int assign[])
{
	double min_dist = 10000000;
	double x=0, y=0, temp_dist=0;
	int k_min_index = 0;

	for(int i = 0; i < (numOfElements/num_of_processes) + 1; i++)
	{
		for(int j = 0; j < numOfClusters; j++)
		{
			x = abs(recv_x[i] - k_x[j]);
			y = abs(recv_y[i] - k_y[j]);
			temp_dist = sqrt((x*x) + (y*y));

			// new minimum distance found
			if(temp_dist < min_dist)
			{
				min_dist = temp_dist;
				k_min_index = j;
			}
		}

		// update the cluster assignment of this data points
		assign[i] = k_min_index;
	}

}
void assign(double k_x[], double k_y[], double recv_x[], double recv_y[], int assign[]){
	double min_dist = 10000000;
	double x=0, y=0, temp_dist=0;
	int k_min_index = 0;

	for(int i = 0; i < (numOfElements); i++)
	{
		for(int j = 0; j < numOfClusters; j++)
		{
			x = abs(recv_x[i] - k_x[j]);
			y = abs(recv_y[i] - k_y[j]);
			temp_dist = sqrt((x*x) + (y*y));

			// new minimum distance found
			if(temp_dist < min_dist)
			{
				min_dist = temp_dist;
				k_min_index = j;
			}
		}

		// update the cluster assignment of this data points
		assign[i] = k_min_index;
	}

}

/* Recalcuate k-means of each cluster because each data point may have
   been reassigned to a new cluster for each iteration of the algorithm */
void calcKmeans(double k_means_x[], double k_means_y[], double data_x_points[], double data_y_points[], int k_assignment[])
{
	double total_x = 0;
	double total_y = 0;
	int numOfpoints = 0;

	for(int i = 0; i < numOfClusters; i++)
	{
		total_x = 0;
		total_y = 0;
		numOfpoints = 0;

		for(int j = 0; j < numOfElements; j++)
		{
			if(k_assignment[j] == i)
			{
				total_x += data_x_points[j];
				total_y += data_y_points[j];
				numOfpoints++;
			}
		}

		if(numOfpoints != 0)
		{
			k_means_x[i] = total_x / numOfpoints;
			k_means_y[i] = total_y / numOfpoints;
		}
	}

}

void out_file(double *x, double *y, int *k,int numOfElements){

    // char dir[105] = "./output/";
    // strcat(dir,cluster_filename);
	FILE *fout = fopen("output/data1.txt", "w");

    for(int i = 0; i < numOfElements; i++){

        fprintf(fout, "%f %f %d\n",x[i],y[i],k[i]);

    }

    fclose(fout);
    system("gnuplot -p -e \"plot 'output/data1.txt' using 1:2:3 with points palette notitle\"");
    // remove("data1.txt");

}

int main(int argc, char *argv[])
{
	// initialize the MPI environment
	MPI_Init(NULL, NULL);

	// get number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// get rank
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// send buffers
	double *k_means_x = NULL;		// k means corresponding x values
	double *k_means_y = NULL;		// k means corresponding y values
	int *k_assignment = NULL;		// each data point is assigned to a cluster
	double *data_x_points = NULL;
	double *data_y_points = NULL;

	// receive buffer
	double *recv_x = NULL;
	double *recv_y = NULL;
	int *recv_assign = NULL;

	if(world_rank == 0)
	{
		printf("%d",world_size);
		if(argc != 2)
		{
			printf("Please include an argument after the program name to list how many processes.\n");
			printf("e.g. To indicate 4 processes, run: mpirun -n 4 ./kmeans 4\n");
			exit(-1);
		}

		num_of_processes = atoi(argv[1]);
		printf("How many clusters would you like to analyze for? ");
		char buffer[2];
		
		scanf("%s", buffer);
		printf("\n");

		numOfClusters = atoi(buffer);
		printf("Ok %d clusters it is.\n", numOfClusters);


		// broadcast the number of clusters to all nodes
		MPI_Bcast(&numOfClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for arrays
		k_means_x = (double *)malloc(sizeof(double) * numOfClusters);
		k_means_y = (double *)malloc(sizeof(double) * numOfClusters);

		if(k_means_x == NULL || k_means_y == NULL)
		{
			perror("malloc");
			exit(-1);
		}

		printf("Reading input data from file...\n\n");

		FILE* fp = fopen("datasets/dataset-10000.txt", "r");

		if(!fp)
		{
			perror("fopen");
			exit(-1);
		}

		// count number of lines to find out how many elements
		// int c = 0;
		// numOfElements = 0;
		FILE *fin = fopen("./datasets/dataset-100000.txt", "r");
        fscanf(fin, "%d", &numOfElements);

		printf("There are a total number of %d elements in the file.\n", numOfElements);

		// broadcast the number of elements to all nodes
		MPI_Bcast(&numOfElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for an array of data points
		data_x_points = (double *)malloc(sizeof(double) * numOfElements);
		data_y_points = (double *)malloc(sizeof(double) * numOfElements);
		k_assignment = (int *)malloc(sizeof(int) * numOfElements);

		if(data_x_points == NULL || data_y_points == NULL || k_assignment == NULL)
		{
			perror("malloc");
			exit(-1);
		}

		for(int i=0;i<numOfElements;i++){
            int x,y;
            fscanf(fin, "%d %d",&x,&y);
            	data_x_points[i] = x;
		    	data_y_points[i] = y;
                k_assignment[i] = 0;
        }

		// close file pointer
		fclose(fp);

		// randomly select initial k-means
		// time_t t;
		// srand((unsigned) time(&t));
		srand(50);
		int random;
		for(int i = 0; i < numOfClusters; i++) {
			random = rand() % numOfElements;
			k_means_x[i] = (rand()-RAND_MAX/2) % (int) max_range;
			k_means_y[i] = (rand()-RAND_MAX/2) % (int) max_range;
		}

		printf("Running k-means algorithm for %d iterations...\n\n", MAX_ITERATIONS);
		for(int i = 0; i < numOfClusters; i++)
		{
			printf("Initial K-means: (%f, %f)\n", k_means_x[i], k_means_y[i]);
		}

		// allocate memory for receive buffers
		recv_x = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes) + 1));
		recv_y = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes) + 1));
		recv_assign = (int *)malloc(sizeof(int) * ((numOfElements/num_of_processes) + 1));

		if(recv_x == NULL || recv_y == NULL || recv_assign == NULL)
		{
			perror("malloc");
			exit(-1);
		}
	}
	else
	{	// I am a worker node
		printf("Worker node");
		num_of_processes = atoi(argv[1]);

		// receive broadcast of number of clusters
		MPI_Bcast(&numOfClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// receive broadcast of number of elements
		MPI_Bcast(&numOfElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for arrays
		k_means_x = (double *)malloc(sizeof(double) * numOfClusters);
		k_means_y = (double *)malloc(sizeof(double) * numOfClusters);

		if(k_means_x == NULL || k_means_y == NULL)
		{
			perror("malloc");
			exit(-1);
		}

		// allocate memory for receive buffers
		recv_x = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes) + 1));
		recv_y = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes) + 1));
		recv_assign = (int *)malloc(sizeof(int) * ((numOfElements/num_of_processes) + 1));

		if(recv_x == NULL || recv_y == NULL || recv_assign == NULL)
		{
			perror("malloc");
			exit(-1);
		}
	}

	/* Distribute the work among all nodes. The data points itself will stay constant and
	   not change for the duration of the algorithm. */
	MPI_Scatter(data_x_points, (numOfElements/num_of_processes) + 1, MPI_DOUBLE,
		recv_x, (numOfElements/num_of_processes) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(data_y_points, (numOfElements/num_of_processes) + 1, MPI_DOUBLE,
		recv_y, (numOfElements/num_of_processes) + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int count = 0;

	// printf("Starting iterate..\n");
    double t1 = MPI_Wtime(); 
	while(count < MAX_ITERATIONS)
	{
		// broadcast k-means arrays
		MPI_Bcast(k_means_x, numOfClusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(k_means_y, numOfClusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// scatter k-cluster assignments array
		MPI_Scatter(k_assignment, (numOfElements/num_of_processes) + 1, MPI_INT,
			recv_assign, (numOfElements/num_of_processes) + 1, MPI_INT, 0, MPI_COMM_WORLD);

		// assign the data points to a cluster
		assign2Cluster(k_means_x, k_means_y, recv_x, recv_y, recv_assign);

		// gather back k-cluster assignments
		MPI_Gather(recv_assign, (numOfElements/num_of_processes)+1, MPI_INT,
			k_assignment, (numOfElements/num_of_processes)+1, MPI_INT, 0, MPI_COMM_WORLD);

		// let the root process recalculate k means
		if(world_rank == 0)
		{
			calcKmeans(k_means_x, k_means_y, data_x_points, data_y_points, k_assignment);
			// printf("Finished iteration %d\n",count);
		}

		count++;
	}

	double t2 = MPI_Wtime(); 
    double duration = t2- t1;

	if(world_rank == 0)
	{
		printf("--------------------------------------------------\n");
		printf("FINAL RESULTS:\n");
		assign(k_means_x, k_means_y, data_x_points, data_y_points, k_assignment);

		out_file(data_x_points, data_y_points, k_assignment,numOfElements);
		for(int i = 0; i < numOfClusters; i++)
		{
			printf("Cluster #%d: (%f, %f)\n", i, k_means_x[i], k_means_y[i]);
		}
		printf("--------------------------------------------------\n");
		printf("Duration: %f\n",duration);
	}

	// deallocate memory and clean up
	free(k_means_x);
	free(k_means_y);
	free(data_x_points);
	free(data_y_points);
	free(k_assignment);
	free(recv_x);
	free(recv_y);
	free(recv_assign);

	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

}
