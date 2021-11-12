#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "must provide exactly 2 arguments!\n");
		return 1;
	}
	int num_td;
	unsigned long long r = atoll(argv[1]);
	unsigned long long k = atoll(argv[2]);
	unsigned long long sum = 0;
	unsigned long long pixels = 0;
	cpu_set_t cpuset;
	sched_getaffinity(0, sizeof(cpuset), &cpuset);
	num_td = CPU_COUNT(&cpuset);
	omp_lock_t lock; 
	omp_init_lock(&lock);
#pragma omp parallel num_threads(num_td) firstprivate(pixels) shared(r, k)
	{
		int id = omp_get_thread_num();
		// printf("id = %d, pixels = %ld\n", id, pixels);
		// std::cout << "id = " << id << " pixels = " << pixels << "\n";
		for (unsigned long long x = id; x < r; x+=num_td) {
			unsigned long long y = ceil(sqrtl(r*r - x*x));
			pixels += y;
		}
		omp_set_lock(&lock);
		sum += pixels;
		omp_unset_lock(&lock);
	}
	sum %= k;
	printf("%llu\n", (4 * sum) % k);
	omp_destroy_lock(&lock);
}
