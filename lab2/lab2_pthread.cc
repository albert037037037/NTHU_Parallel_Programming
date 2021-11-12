#include <pthread.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

unsigned long long r, k;
unsigned long long sum;
unsigned long long num_threads;
pthread_mutex_t mutex;

void* calculate(void* threadid) {
    int* tid = (int*)threadid;	
	unsigned long long pixels = 0;
	for (unsigned long long x = *tid; x < r; x+=num_threads) {
		unsigned long long y = ceil(sqrtl(r*r - x*x));
		pixels += y;
	}
	pthread_mutex_lock(&mutex);
	sum += pixels;
	pthread_mutex_unlock(&mutex);
    pthread_exit(NULL);
}

int main(int argc, char** argv) {
	if (argc != 3) {
		fprintf(stderr, "must provide exactly 2 arguments!\n");
		return 1;
	}
	r = atoll(argv[1]);
	k = atoll(argv[2]);
	pthread_mutex_init (&mutex, NULL);
	cpu_set_t cpuset;
	sched_getaffinity(0, sizeof(cpuset), &cpuset);
	num_threads = CPU_COUNT(&cpuset);
	pthread_t threads[num_threads];
	sum = 0;
    int rc;
    int ID[num_threads];
    int t;
    for (t = 0; t < num_threads; t++) {
        ID[t] = t;
        rc = pthread_create(&threads[t], NULL, calculate, (void*)&ID[t]);
    }

	for (int i=0; i<num_threads; i++) {
		pthread_join(threads[i], NULL);
	}
	
	sum %= k;
	pthread_mutex_destroy(&mutex);
	printf("%llu\n", (4 * sum) % k);
}
