#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv)
{
	int rank, size;
	double starttime, endtime;
	MPI_Init(&argc, &argv);
	starttime = MPI_Wtime();
	if (argc != 3)
	{
		fprintf(stderr, "must provide exactly 2 arguments!\n");
		return 1;
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	unsigned long long r = atoll(argv[1]);
	unsigned long long k = atoll(argv[2]);
	unsigned long long pixels = 0;
	unsigned long long sum;
	// unsigned long long count = r / size;
	// if(rank != size-1) {
	// 	for (unsigned long long x = rank * count; x < rank*count+count; x += 1)
	// 	{
	// 		unsigned long long y = ceil(sqrtl((r * r - x * x)));
	// 		pixels += y;
	// 		pixels %= k;
	// 	}
	// }
	// else {
	// 	for (unsigned long long x = rank * count; x < r; x += 1)
	// 	{
	// 		unsigned long long y = ceil(sqrtl((r * r - x * x)));
	// 		pixels += y;
	// 		pixels %= k;
	// 	}
	// }
 	for (unsigned long long x = rank; x < r; x += size)
	{
		unsigned long long y = ceil(sqrtl((r * r - x * x)));
		pixels += y;
		pixels %= k;
	}

	MPI_Reduce(&pixels, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, size - 1, MPI_COMM_WORLD);
	if (rank == size - 1)
		printf("%llu\n", (4 * sum) % k);
	endtime = MPI_Wtime();
	printf("That took %f seconds\n", endtime - starttime);
	MPI_Finalize();
}
