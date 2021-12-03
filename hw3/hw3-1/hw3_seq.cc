#include <bits/stdc++.h>
#include <climits>
#include <pthread.h>

#define MAX_V 6000 // In cpu version, max vertices are 6000
int V, E;

void floydWarshall(int **graph)
{
	int dist[V][V], i, j, k;
	for (i = 0; i < V; i++)
		for (j = 0; j < V; j++)
			dist[i][j] = graph[i][j];

	for (k = 0; k < V; k++) {
		for (i = 0; i < V; i++) {
			for (j = 0; j < V; j++) {
				if (dist[i][j] > (dist[i][k] + dist[k][j]) && (dist[k][j] != INT_MAX && dist[i][k] != INT_MAX)) {
					dist[i][j] = dist[i][k] + dist[k][j];
                }
			}
		}
	}
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
            else graph[i][j] = INT_MAX;
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

int main(int argc, char** argv)
{
    assert(argc == 2);
    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];
    int **graph = readInputFile(inputFileName);
	floydWarshall(graph);
	return 0;
}


