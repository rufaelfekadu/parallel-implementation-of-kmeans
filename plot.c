#include <stdlib.h>
#include <stdio.h>
int main(){
        system("gnuplot -p -e \"plot 'cluster_output_dataset1.txt' using 1:2:3 with points palette notitle\"");
return 0;
}