#include <stdlib.h>
#include <stdio.h>
int main(){
        system("gnuplot -p -e \"splot 'cluster_output_dataset5.txt' using 1:2:3:4 with points palette notitle\"");
return 0;
}