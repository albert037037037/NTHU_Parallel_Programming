#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>
#include <mpi.h>

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "must provide exactly 2 arguments!\n");
		return 1;
	}
	MPI_Init(&argc, &argv);
    int rank, size, omp_threads, omp_thread;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	int num_td;
	unsigned long long r = atoll(argv[1]);
	unsigned long long k = atoll(argv[2]);
	unsigned long long local_sum = 0;
	unsigned long long total_sum = 0;
	unsigned long long pixels = 0;

	cpu_set_t cpuset;
	sched_getaffinity(0, sizeof(cpuset), &cpuset);
	num_td = CPU_COUNT(&cpuset);
	omp_lock_t lock; 
	omp_init_lock(&lock);
#pragma omp parallel num_threads(num_td) firstprivate(pixels) shared(local_sum, r, k, rank, size, num_td)
	{
		int id = omp_get_thread_num();
		for (unsigned long long x = rank*num_td+id; x < r; x+=(size*num_td)) {
			unsigned long long y = ceil(sqrtl(r*r - x*x));
			pixels += y;
		}
		omp_set_lock(&lock);
		local_sum += pixels;
		omp_unset_lock(&lock);
	}
	MPI_Reduce(&local_sum, &total_sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	total_sum %= k;
	if(rank == 0 )printf("%llu\n", (4 * total_sum) % k);
	omp_destroy_lock(&lock);
}
