
serial:
	gcc -o output/compiled/serial serial.c -lm -fopenmp

openmp:
	gcc -o output/compiled/openmp parallel_openmp.c -fopenmp -lm

pthread:
	gcc -o output/compiled/pthread parallel_pthread.c -fopenmp -pthread -lm

mpi:
	mpicxx -o output/compiled/mpi parallel_mpi.c -lm

mpi_openmp:
	mpicxx -o output/compiled/mpi_openmp parallel_mpi_with_openmp.c -fopenmp -lm

plot:
	gcc -o plot plot.c

clean:
	$(RM) output/compiled/*;

clean_data:
	$(RM) output/to_plot/*;
	$(RM) output/plot/*;


