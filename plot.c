#include <stdlib.h>
#include <stdio.h>

void plot_vs_cluster(){

    system("gnuplot -e \"set xlabel 'number of clusters'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/plots/time_vs_num_clusters.png'; plot 'output/to_plot/serial_vs_cluster.txt' using 1:2 title 'Serial' with lines, 'output/to_plot/openmp_vs_cluster.txt' using 1:2 title 'OpenMP' with lines, 'output/to_plot/mpi_vs_cluster.txt' using 1:2 title 'MPI' with lines, 'output/to_plot/pthread_vs_cluster.txt' using 1:2 title 'pthread' with lines, 'output/to_plot/mpi_openmp_vs_cluster.txt' using 1:2 title 'MPI_OpenMP' with lines lw 2 \"");

}
void plot_vs_thread(){

    system("gnuplot -e \"set xlabel 'number of threads'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/plots/time_vs_num_threads.png'; plot 'output/to_plot/openmp_vs_thread.txt' using 1:2 title 'OpenMP' with lines, 'output/to_plot/pthread_vs_thread.txt' using 1:2 title 'pthread' with lines, 'output/to_plot/mpi_openmp_vs_thread.txt' using 1:2 title 'MPI_OpenMP' with lines lw 2\"");

}
void plot_vs_process(){

    system("gnuplot -e \"set xlabel 'number of threads'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/plots/time_vs_num_procs.png'; plot 'output/to_plot/mpi_vs_process.txt' using 1:2 title 'OpenMP' with lines, 'output/to_plot/mpi_openmp_vs_process.txt' using 1:2 title 'pthread' with lines lw 2\"");

}
void plot_vs_point(){

    system("gnuplot -e \"set xlabel 'number of points'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/plots/time_vs_num_points.png'; plot 'output/to_plot/serial_vs_point.txt' using 1:2 title 'Serial' with lines, 'output/to_plot/openmp_vs_point.txt' using 1:2 title 'OpenMP' with lines, 'output/to_plot/mpi_vs_point.txt' using 1:2 title 'MPI' with lines, 'output/to_plot/pthread_vs_point.txt' using 1:2 title 'pthread' with lines, 'output/to_plot/mpi_openmp_vs_point.txt' using 1:2 title 'MPI_OpenMP' with lines lw 2 \"");

}
void plot_data(){

    system("gnuplot -p -e \"plot 'datasets/dataset-10000.txt' using 1:2 with points notitle\"");

}
int main(){
        plot_vs_point();
        plot_vs_cluster();
        plot_vs_thread();
        plot_vs_process();
        // plot_data();
        // system("gnuplot -p -e \"splot 'cluster_output_dataset5.txt' using 1:2:3:4 with points palette notitle\"");
return 0;
}