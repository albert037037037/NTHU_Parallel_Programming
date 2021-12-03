#include <bits/stdc++.h>
#include <climits>
#include <omp.h>
#include <pthread.h>

pthread_barrier_t barrier;
pthread_mutex_t mutexUse;

#define INF 1073741823 // In cpu version, max vertices are 6000
int V, E;

typedef struct _data {
    int numThreads;
    int** graph;
    int* nowUse;
} myData;

void* floydWarshall(void* void_data)
{
    myData* data = (myData*) void_data;
    pthread_mutex_lock(&mutexUse);
        int useI = *(data->nowUse);
        *(data->nowUse) += 1;
    pthread_mutex_unlock(&mutexUse);
	for (int k = 0; k < V; k++) {
		for (int i = useI; i < V; i+=data->numThreads) {
			for (int j = 0; j < V; j++) {
				if (data->graph[i][j] > (data->graph[i][k] + data->graph[k][j]) && (data->graph[k][j] != INF && data->graph[i][k] != INF)) {
					data->graph[i][j] = data->graph[i][k] + data->graph[k][j];
                }
			}
		}
        pthread_barrier_wait(&barrier);
	}
    pthread_exit(data);
}

int** readInputFile(const char* fileName) {
    FILE* fp = fopen(fileName, "r");
    fread(&V, sizeof(int), 1, fp);
    fread(&E, sizeof(int), 1, fp);

    int** graph = (int**) malloc(sizeof(int*)*V);
    for(int i=0; i<V; i++) {
        graph[i] = (int*) malloc(sizeof(int)*V);
        for(int j=0; j<V; j++) {
            if(i == j) graph[i][j] = 0;
            else graph[i][j] = INF;
        }
    }
    
    int srcV, destV, weight;
    for(int i=0; i<E; i++) {
        fread(&srcV, sizeof(int), 1, fp);
        fread(&destV, sizeof(int), 1, fp);
        fread(&weight, sizeof(int), 1, fp);
        graph[srcV][destV] = weight;
    }
    fclose(fp);
    return graph;
}

void writeOutputFile(const char* fileName, myData* data) {
    FILE* fp = fopen(fileName, "wb");
    for(int i=0; i<V; i++) {
        fwrite(data->graph[i], sizeof(int), V, fp);
    }
    fclose(fp);
}

int main(int argc, char** argv)
{
    assert(argc == 3);
    // check cpu
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    int numThreads = CPU_COUNT(&cpu_set);
    pthread_t threads[numThreads];
    int ID[numThreads];

    // read input file
    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];
    int **graph = readInputFile(inputFileName);

    // initial pthread data
    myData* data = (myData*) malloc(sizeof(myData));
    data->graph = graph;
    data->numThreads = numThreads;
    data->nowUse = (int*) malloc(sizeof(int));
    *(data->nowUse) = 0;

    // initial lock and barrier
    pthread_mutex_init (&mutexUse, NULL);
    pthread_barrier_init(&barrier, NULL, numThreads);

    // pthread function
    for (int t = 0; t < numThreads; t++) {
        ID[t] = t;
        pthread_create(&threads[t], NULL, floydWarshall, (void*)data);
    }

    // pthread join
    for (int i=0; i<numThreads; i++) {
		pthread_join(threads[i], (void**)&data);
	}

    // output file
    writeOutputFile(outputFileName, data);

    // free memory
    free(data->nowUse);
    free(graph);
    free(data);
    pthread_mutex_destroy(&mutexUse);
	return 0;
}


