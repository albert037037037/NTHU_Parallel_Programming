#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int INF = ((1 << 30) - 1);
void input(char* inFileName);
void output(char* outFileName);

void block_FW();
int ceil(int a, int b);
void cal(int B, int Round, int block_start_x, int block_start_y, int block_width, int block_height);

int n, m;
int padN;
int* Dist;

int main(int argc, char* argv[]) {
    input(argv[1]);
    block_FW();
    output(argv[2]);
    return 0;
}

void input(char* infile) {
    FILE* file = fopen(infile, "rb");
    fread(&n, sizeof(int), 1, file);
    fread(&m, sizeof(int), 1, file);
    printf("vertice = %d\n", n);
    int pad = n % 32;
    padN = n + (32-pad);
    Dist = (int*) malloc(sizeof(int)*padN*padN);
    for (int i = 0; i < padN; ++i) {
        for (int j = 0; j < padN; ++j) {
            if (i == j) {
                Dist[i*padN+j] = 0;
            } else {
                Dist[i*padN+j] = INF;
            }
        }
    }

    int pair[3];
    for (int i = 0; i < m; ++i) {
        fread(pair, sizeof(int), 3, file);
        Dist[pair[0]*padN+pair[1]]= pair[2];
    }
    fclose(file);
}

void output(char* outFileName) {
    FILE* outfile = fopen(outFileName, "w");
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (Dist[i*padN+j] >= INF) Dist[i*padN+j] = INF;
        }
        fwrite(&Dist[i*padN], sizeof(int), n, outfile);
    }
    fclose(outfile);
}

__global__ void blockedASAPPhase1(int round, int* graph, int V, int blockSize) {
    __shared__ int shareDistInBlock[32*32];
    int getBlockSize = blockSize;
    int i = blockSize * round + threadIdx.y;
    int j = blockSize * round + threadIdx.x;
    shareDistInBlock[threadIdx.y * getBlockSize + threadIdx.x] = graph[i*V+j];
    __syncthreads();
    for(int k = 0; k<getBlockSize; k++) {
        shareDistInBlock[threadIdx.y * getBlockSize + threadIdx.x] = min(shareDistInBlock[threadIdx.y * getBlockSize + threadIdx.x], shareDistInBlock[threadIdx.y*getBlockSize + k] + shareDistInBlock[k*getBlockSize + threadIdx.x]);
        // if(shareDistInBlock[threadIdx.x * getBlockSize + threadIdx.y] > shareDistInBlock[threadIdx.x*getBlockSize + k] + shareDistInBlock[k*getBlockSize + threadIdx.y]){
        //     shareDistInBlock[threadIdx.x * getBlockSize + threadIdx.y] = shareDistInBlock[threadIdx.x*getBlockSize + k] + shareDistInBlock[k*getBlockSize + threadIdx.y];
        // }
        __syncthreads();
    }
    graph[i*V+j] = shareDistInBlock[threadIdx.y*getBlockSize+threadIdx.x];
}

// __global__ void blockedASAPPhase2Row(int round, int* graph, int V, int blockSize) {
//     if(blockIdx.x != round) {
//         __shared__ int shareDistInBlockPivot[32*32];
//         __shared__ int shareDistInBlockOther[32*32];
//         int getBlockSize = blockSize;
//         int pivotI = getBlockSize * round + threadIdx.x;
//         int pivotJ = getBlockSize * round + threadIdx.y;
//         int targetI = pivotI;
//         int targetJ = blockIdx.x * getBlockSize + threadIdx.y;

//         if(pivotI < V && pivotJ < V) {
//             shareDistInBlockPivot[threadIdx.x * getBlockSize + threadIdx.y] = graph[pivotI*V+pivotJ];
//         }
//         else {
//             shareDistInBlockPivot[threadIdx.x * getBlockSize + threadIdx.y] = INF;
//         }

//         if(targetI < V && targetJ < V) {
//             shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] = graph[targetI*V+targetJ];
//         }
//         else {
//             shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] = INF;
//         }
//         __syncthreads();
//         if(targetI < V && targetJ < V) {
//             for(int k=0; k<getBlockSize; k++) {
//                 if(shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] > shareDistInBlockPivot[threadIdx.x*getBlockSize + k] + shareDistInBlockOther[k*getBlockSize + threadIdx.y]) {
//                     shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] = shareDistInBlockPivot[threadIdx.x*getBlockSize + k] + shareDistInBlockOther[k*getBlockSize + threadIdx.y];
//                 }
//             }
//             graph[targetI*V+targetJ] = shareDistInBlockOther[threadIdx.x*getBlockSize+threadIdx.y];
//         }
//     }
// }

// __global__ void blockedASAPPhase2Col(int round, int* graph, int V, int blockSize) {
//     if(blockIdx.x != round) {
//         __shared__ int shareDistInBlockPivot[32*32];
//         __shared__ int shareDistInBlockOther[32*32];
//         int getBlockSize = blockSize;
//         int pivotI = getBlockSize * round + threadIdx.x;
//         int pivotJ = getBlockSize * round + threadIdx.y;
//         int targetI = blockIdx.x * getBlockSize + threadIdx.x;
//         int targetJ = pivotJ;

//         if(pivotI < V && pivotJ < V) {
//             shareDistInBlockPivot[threadIdx.x * getBlockSize + threadIdx.y] = graph[pivotI*V+pivotJ];
//         }
//         else {
//             shareDistInBlockPivot[threadIdx.x * getBlockSize + threadIdx.y] = INF;
//         }

//         if(targetI < V && targetJ < V) {
//             shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] = graph[targetI*V+targetJ];
//         }
//         else {
//             shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] = INF;
//         }
//         __syncthreads();
//         if(targetI < V && targetJ < V) {
//             for(int k=0; k<getBlockSize; k++) {
//                 if(shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] > shareDistInBlockOther[threadIdx.x*getBlockSize + k] + shareDistInBlockPivot[k*getBlockSize + threadIdx.y]) {
//                     shareDistInBlockOther[threadIdx.x * getBlockSize + threadIdx.y] = shareDistInBlockOther[threadIdx.x*getBlockSize + k] + shareDistInBlockPivot[k*getBlockSize + threadIdx.y];
//                 }
//             }
//             graph[targetI*V+targetJ] = shareDistInBlockOther[threadIdx.x*getBlockSize+threadIdx.y];
//         }
//     }
// }

__global__ void blockedASAPPhase2(int round, int* graph, int V, int blockSize) {
    if(blockIdx.x != round) {
        __shared__ int shareDistInBlockPivot[32*32];
        __shared__ int shareDistInBlockOtherRow[32*32];
        __shared__ int shareDistInBlockOtherCol[32*32];
        int getBlockSize = blockSize;
        int pivotI = getBlockSize * round + threadIdx.y;
        int pivotJ = getBlockSize * round + threadIdx.x;
        int targetIRow = pivotI;
        int targetJRow = blockIdx.x * getBlockSize + threadIdx.x;
        int targetICol = blockIdx.x * getBlockSize + threadIdx.y;
        int targetJCol = pivotJ;
        shareDistInBlockPivot[threadIdx.y * getBlockSize + threadIdx.x] = graph[pivotI*V+pivotJ];
        shareDistInBlockOtherRow[threadIdx.y * getBlockSize + threadIdx.x] = graph[targetIRow*V+targetJRow];
        shareDistInBlockOtherCol[threadIdx.y * getBlockSize + threadIdx.x] = graph[targetICol*V+targetJCol];
        __syncthreads();
        for(int k=0; k<getBlockSize; k++) {
            shareDistInBlockOtherRow[threadIdx.y * getBlockSize + threadIdx.x] = min(shareDistInBlockOtherRow[threadIdx.y * getBlockSize + threadIdx.x], shareDistInBlockPivot[threadIdx.y*getBlockSize + k] + shareDistInBlockOtherRow[k*getBlockSize + threadIdx.x]);
            shareDistInBlockOtherCol[threadIdx.y * getBlockSize + threadIdx.x] = min(shareDistInBlockOtherCol[threadIdx.y * getBlockSize + threadIdx.x], shareDistInBlockOtherCol[threadIdx.y*getBlockSize + k] + shareDistInBlockPivot[k*getBlockSize + threadIdx.x]);
        }
        graph[targetIRow*V+targetJRow] = shareDistInBlockOtherRow[threadIdx.y*getBlockSize+threadIdx.x];
        graph[targetICol*V+targetJCol] = shareDistInBlockOtherCol[threadIdx.y*getBlockSize+threadIdx.x];
    }
}

__global__ void blockedASAPPhase3(int round, int* graph, int V, int blockSize) {
    if(blockIdx.x != round && blockIdx.y != round) {
        __shared__ int shareDistInBlockPivotRow[32*32];
        __shared__ int shareDistInBlockPivotCol[32*32];
        int getBlockSize = blockSize;
        int pivotRowI = getBlockSize * round + threadIdx.y;
        int pivotRowJ = blockIdx.y * blockDim.x + threadIdx.x;
        int pivotColI = blockIdx.x * blockDim.y + threadIdx.y;
        int pivotColJ = getBlockSize * round + threadIdx.x;
        int pointUpdatedByThisThread = graph[ pivotColI*V + pivotRowJ];
        shareDistInBlockPivotRow[threadIdx.y * getBlockSize + threadIdx.x] = graph[pivotRowI*V+pivotRowJ];
        shareDistInBlockPivotCol[threadIdx.y * getBlockSize + threadIdx.x] = graph[pivotColI*V+pivotColJ];
        __syncthreads();

        for(int k=0; k<getBlockSize; k++) {
            pointUpdatedByThisThread = min(pointUpdatedByThisThread, shareDistInBlockPivotCol[threadIdx.y*getBlockSize + k] + shareDistInBlockPivotRow[k*getBlockSize + threadIdx.x]);
            // if(pointUpdatedByThisThread > shareDistInBlockPivotCol[threadIdx.x*getBlockSize + k] + shareDistInBlockPivotRow[k*getBlockSize + threadIdx.y]) {
            //     pointUpdatedByThisThread = shareDistInBlockPivotCol[threadIdx.x*getBlockSize + k] + shareDistInBlockPivotRow[k*getBlockSize + threadIdx.y];
            // }
        }
        graph[pivotColI*V + pivotRowJ] = pointUpdatedByThisThread;
    }
}

int ceil(int a, int b) { return (a + b - 1) / b; }

void block_FW() {
    // malloc space on gpu
    int * deviceGraph;
    deviceGraph = (int *)malloc(padN*padN*sizeof(int));
    cudaHostRegister(Dist, padN*padN*sizeof(int), cudaHostRegisterDefault);
    cudaMalloc(&deviceGraph, padN*padN*sizeof(int));
    cudaMemcpy(deviceGraph, Dist,  padN*padN*sizeof(int), cudaMemcpyHostToDevice);

    // initial parameter
    int round = ceil(n, 32);
    dim3 gridSize(1, 1);
    dim3 gridSize2(round, 1);
    dim3 gridSize3(round, round);
    dim3 blockSize(32, 32);
    
    // for(int i=0; i<n; i++) {
    //     for(int j=0; j<n; j++) {
    //         printf("%d ", Dist[i*n+j]);
    //     }
    //     printf("\n");
    // }

    for (int r = 0; r < round; ++r) {
        printf("%d %d\n", r, round);
        fflush(stdout);
        /* Phase 1*/
        // calculate diagonal axis block (upper left to down right)
        blockedASAPPhase1<<<gridSize, blockSize>>> (r, deviceGraph, padN, 32);

        /* Phase 2*/
        // calculate horizontal axis block
        // blockedASAPPhase2Row<<<gridSize2, blockSize>>> (r, deviceGraph, n, 32);
        // calculate vertical axis block
        // blockedASAPPhase2Col<<<gridSize2, blockSize>>> (r, deviceGraph, n, 32);
        blockedASAPPhase2<<<gridSize2, blockSize>>> (r, deviceGraph, padN, 32);
        /* Phase 3*/
        blockedASAPPhase3<<<gridSize3, blockSize>>> (r, deviceGraph, padN, 32);
    }
    cudaMemcpy(Dist, deviceGraph,  padN*padN*sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(deviceGraph);
    // for(int i=0; i<n; i++) {
    //     for(int j=0; j<n; j++) {
    //         printf("%d ", Dist[i*n+j]);
    //     }
    //     printf("\n");
    // }
    cudaFree(deviceGraph);
}
