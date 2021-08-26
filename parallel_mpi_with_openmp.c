#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <assert.h>
#include <float.h>

#define MATCH(s) (!strcmp(argv[ac], (s)))

int num_elements,max_iterations=1000;
int max_range = 800;
// Creates an array of random floats. Each number has a value from 0 - 1
float* init_points(const char *dataset_filename,int site_per_proc,int nprocs) {

  float *all_points = (float *)malloc(sizeof(float) * num_elements * 2);
  int * num;
  char dir[105]="./datasets/";
  strcat(dir,dataset_filename); 
	FILE *fin = fopen(dir, "r");
  fscanf(fin, "%d", &num_elements);

  for(int i=0;i<num_elements*2;i++){

      fscanf(fin, "%f",&all_points[i]);
  }

  fclose(fin);

  return all_points;
}

// Distance**2 between d-vectors pointed to by v1, v2.
float distance2(const float *v1, const float *v2, const int d) {
  float dist = 0.0;
  for (int i=0; i<d; i++) {
    float diff = v1[i] - v2[i];
    dist += diff * diff;
  }
  return dist;
}

// Assign a site to the correct cluster by computing its distances to
// each cluster centroid.
int assign_site(const float* site, float* centroids, const int k, const int d, const int numthreads) {
  int best_cluster = 0;
  float best_dist = distance2(site, centroids, d);
  float* centroid = centroids + d;
  int c ;
 
    #pragma omp for  schedule(static)
    for (c = 1; c < k; c++) {
        
        float dist = distance2(site, centroid, d);
        if (dist < best_dist) {
            best_cluster = c;
            best_dist = dist;
        }
        centroid += d;
    }
  
  return best_cluster;
}


// Add a site (vector) into a sum of sites (vector).
void add_site(const float * site, float * sum, const int d) {
  for (int i=0; i<d; i++) {
  sum[i] += site[i];
  }
}

// Print the centroids one per line.
void print_centroids(float * centroids, const int k, const int d) {
  float *p = centroids;
  printf("Centroids:\n");
  for (int i = 0; i<k; i++) {
    for (int j = 0; j<d; j++, p++) {
      printf("%f ", *p);
    }
    printf("\n");
  }
}

void mpi_openmp_vs_cluster(double dur,int num_cluster){
    char dir[105]="./output/to_plot/mpi_openmp_vs_cluster.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",num_cluster, dur);

    fclose(fout);
}

void mpi_openmp_vs_point(double dur,int size){
    char dir[105]="./output/to_plot/mpi_openmp_vs_point.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",size, dur);

    fclose(fout);
}

void mpi_openmp_vs_thread(double dur,int numthread){
    char dir[105]="./output/to_plot/mpi_openmp_vs_thread.txt";
	FILE *fout = fopen(dir, "a");
	fprintf(fout, "%d %f\n",numthread, dur);

    fclose(fout);
}

int main(int argc, char** argv) {

    srand(31359);
    int ac = 1;
    int disable_display = 0;
    int seedVal = 100;
    int nx = 1;
    int num_cluster = 4;
    int num_threads = 1;

    for(ac=1;ac<argc;ac++)
    {
        if(MATCH("-n")) {nx = atoi(argv[++ac]);}
        else if(MATCH("-i")) {max_iterations = atoi(argv[++ac]);}
        else if(MATCH("-t"))  {num_threads = atof(argv[++ac]);}
        else if(MATCH("-c"))  {num_cluster = atof(argv[++ac]);}
        // else if(MATCH("-s"))  {seedVal = atof(argv[++ac]);}
        else if(MATCH("-d"))  {disable_display = 1;}
        else {
            printf("Usage: %s [-n < meshpoints>] [-i <iterations>] [-s seed] [-p prob] [-t numthreads] [-step] [-g <game #>] [-d]\n",argv[0]);
            return(-1);
        }
    }

    char dataset_filename[105]="";
    switch (nx)
    {
    case 1:
        strcat(dataset_filename, "dataset-10000.txt");
        num_elements = 10000;
        break;
    case 2:
        strcat(dataset_filename, "dataset-50000.txt");
        num_elements = 50000;
        break;
    case 3:
        strcat(dataset_filename, "dataset-100000.txt");
        num_elements = 100000;
        break;
    case 4:
        strcat(dataset_filename, "dataset-200000.txt");
        num_elements = 200000;
        break;
    case 5:
        strcat(dataset_filename, "dataset-400000.txt");
        num_elements = 400000;
        break;
    case 6:
        strcat(dataset_filename, "dataset-500000.txt");
        num_elements = 500000;
        break;
    case 7:
        strcat(dataset_filename, "dataset-600000.txt");
        num_elements = 600000;
        break;
    case 8:
        strcat(dataset_filename, "dataset-800000.txt");
        num_elements = 800000;
        break;
    case 9:
        strcat(dataset_filename, "dataset-1000000.txt");
        num_elements = 1000000;
        break;
    default:
        strcat(dataset_filename, "dataset-10000.txt");
        num_elements = 10000;
        break;
    }

  // Initial MPI and find process rank and number of processes.
  MPI_Init(NULL, NULL);
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int sites_per_proc = num_elements/nprocs ;
  //
  // Data structures in all processes.
  //
  // The sites assigned to this process.
  float* sites;  

  sites = (float *)malloc(sites_per_proc * 2 * sizeof(float));
  // The sum of sites assigned to each cluster by this process.
  // k vectors of d elements.
  float* sums;

  sums = (float *)malloc(num_cluster * 2 * sizeof(float));
  // The number of sites assigned to each cluster by this process. k integers.
  int* counts;
  
  counts =(int *) malloc(num_cluster * sizeof(int));
  // The current centroids against which sites are being compared.
  // These are shipped to the process by the root process.
  float* centroids;
  
  centroids = (float *)malloc(num_cluster * 2 * sizeof(float));

  // The cluster assignments for each site.
  int* labels;
  
  labels = (int *)malloc(sites_per_proc * sizeof(int));

  // All the sites for all the processes.
  // site_per_proc * nprocs vectors of d floats.
  float* all_sites = NULL;
  // Sum of sites assigned to each cluster by all processes.
  float* grand_sums = NULL;
  // Number of sites assigned to each cluster by all processes.
  int* grand_counts = NULL;
  // Result of program: a cluster label for each site.
  int* all_labels;
  if (rank == 0) {
    float *x,*y;
    all_sites = init_points(dataset_filename,sites_per_proc,nprocs);
    // Take the first k sites as the initial cluster centroids.
    for (int i = 0; i < num_cluster * 2 -1; i++) {
      centroids[i] = all_sites[i];
      // centroids[i] = (rand()-RAND_MAX/2) % (int) max_range; 
    }
    print_centroids(centroids, num_cluster, 2);

    grand_sums = (float *)malloc(num_cluster * 2 * sizeof(float));

    grand_counts = (int *)malloc(num_cluster * sizeof(int));

    all_labels = (int *)malloc(nprocs * sites_per_proc * sizeof(int));
  }

  // Root sends each process its share of sites.
  MPI_Scatter(all_sites,2*sites_per_proc, MPI_FLOAT, sites, 2*sites_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);


  float norm = 1.0;  // Will tell us if centroids have moved.

  float t1 = MPI_Wtime();
  int itr = 1;
  while (norm > 0.00001) { 

    // Broadcast the current cluster centroids to all processes.
    MPI_Bcast(centroids, num_cluster*2, MPI_FLOAT,0, MPI_COMM_WORLD);

    // Each process reinitializes its cluster accumulators.
    for (int i = 0; i < num_cluster*2; i++) sums[i] = 0.0;
    for (int i = 0; i < num_cluster; i++) counts[i] = 0;

    // Find the closest centroid to each site and assign to cluster.
    float* site = sites;
    for (int i = 0; i < sites_per_proc; i++, site += 2) {
      int cluster = assign_site(site, centroids, num_cluster, 2,num_threads);
      // Record the assignment of the site to the cluster.
      counts[cluster]++;
      add_site(site, &sums[cluster*2], 2);
    }

    // Gather and sum at root all cluster sums for individual processes.
    MPI_Reduce(sums, grand_sums, num_cluster * 2, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(counts, grand_counts, num_cluster, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      // Root process computes new centroids by dividing sums per cluster
      // by count per cluster.
      for (int i = 0; i<num_cluster; i++) {
      for (int j = 0; j<2; j++) {
      int dij = 2*i + j;
      grand_sums[dij] /= grand_counts[i];
      }
    }
    // Have the centroids changed much?
    norm = distance2(grand_sums, centroids, 2*num_cluster);
    printf("norm: %f\n",norm);
    // Copy new centroids from grand_sums into centroids.
    for (int i=0; i<num_cluster*2; i++) {
      centroids[i] = grand_sums[i];
    }
    print_centroids(centroids,num_cluster,2);
    }
    // Broadcast the norm.  All processes will use this in the loop test.
    MPI_Bcast(&norm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    itr++;
  }

  // Now centroids are fixed, so compute a final label for each site.
  float* site = sites;
  for (int i = 0; i < sites_per_proc; i++, site += 2) {
    labels[i] = assign_site(site, centroids, num_cluster, 2,num_threads);
  }

  // Gather all labels into root process.
  MPI_Gather(labels, sites_per_proc, MPI_INT, all_labels, sites_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

  float t2 = MPI_Wtime();
  // Root can print out all sites and labels.
  if ((rank == 0) && 1) {
    float time_taken = t2-t1;
    printf("number of iteration: %d\n",itr); 
    printf("time taken per iteration: %f\n",time_taken/itr);

    mpi_openmp_vs_cluster(time_taken,num_cluster);
    mpi_openmp_vs_point(time_taken,num_elements);
    mpi_openmp_vs_thread(time_taken,num_threads);
    if(!disable_display){

      FILE *fout = fopen("output/data.txt", "w");
      float* site = all_sites; 
      for (int i = 0; i < nprocs * sites_per_proc; i++, site += 2) {

        for (int j = 0; j < 2; j++) {
          fprintf(fout,"%f ", site[j]);
        }

        fprintf(fout, "%4d\n",all_labels[i]);
      }
        fclose(fout);
        system("gnuplot -p -e \"plot 'output/data.txt' using 1:2:3 with points palette notitle\"");
        remove("data.txt");
    }

  } 
    

  MPI_Finalize();
  return 0;
}
