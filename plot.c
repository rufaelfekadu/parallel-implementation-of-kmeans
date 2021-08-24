#include <stdlib.h>
#include <stdio.h>

void plot_vs_cluster(){

    system("gnuplot -e \"set xlabel 'number of clusters'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/time_vs_num_clusters.png'; plot 'output/result.txt' using 1:2 title 'Serial' with lines, 'output/result.txt' using 1:3 title 'OpenMP' with lines lw 2 \"");

}
void plot_vs_thread(){

    system("gnuplot -e \"set xlabel 'number of threads'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/time_vs_num_threads.png'; plot 'output/time_vs_num_threads.txt' using 1:2 title 'OpenMP' with lines \"");

}
void plot_vs_point(){

    system("gnuplot -e \"set xlabel 'number of points'; set ylabel 'Execution time'; set terminal png size 800,600; set output 'output/time_vs_num_points.png'; plot 'output/time_vs_num_points.txt' using 1:2 title 'OpenMP' with lines, 'output/time_vs_num_points.txt' using 1:3 title 'Serial' with lines lw 2 \"");

}
void plot_data(){

    system("gnuplot -p -e \"plot 'datasets/dataset-10000.txt' using 1:2 with points notitle\"");

}
int main(){
        // plot_vs_point();
        // plot_vs_cluster();
        // plot_vs_thread();
        plot_data();
        // system("gnuplot -p -e \"splot 'cluster_output_dataset5.txt' using 1:2:3:4 with points palette notitle\"");
return 0;
}