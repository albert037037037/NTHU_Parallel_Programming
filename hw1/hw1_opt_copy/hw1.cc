#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <mpi.h>

#define EvenPhase 998
#define OddPhase 999
float *tmp;

inline void radix_sort(float *ar, int len)
{
    unsigned  *a=(unsigned *)ar;
    unsigned *BUF = (unsigned *)malloc(sizeof(unsigned ) * (len+1)),*b=BUF+1;
    *(float *)BUF=1.0;//guard
    int sum[256] = {0}, sum1[256] = {0}, sum2[256] = {0}, sum3[256] = {0};
    for (int i = 0; i < len; i++)
    {
        ++sum[a[i] & 255];
        ++sum1[(a[i] >> 8) & 255];
        ++sum2[(a[i] >> 16) & 255];
        ++sum3[a[i] >> 24];
    }
    for (int q = 1; q <= 255; ++q)
    {
        sum[q] += sum[q - 1];
        sum1[q] += sum1[q - 1];
        sum2[q] += sum2[q - 1];
        sum3[q] += sum3[q - 1];
    }
    for (int q = len - 1; q >= 0; --q)
        b[--sum[a[q] & 255]] = a[q];
    for (int q = len - 1; q >= 0; --q)
        a[--sum1[(b[q] >> 8) & 255]] = b[q];
    for (int q = len - 1; q >= 0; --q)
        b[--sum2[(a[q] >> 16) & 255]] = a[q];
    for (int q = len - 1; q >= 0; --q)
        a[--sum3[b[q] >> 24]] = b[q];
    memcpy(b,a,sizeof(unsigned)*len);
    int beg=len-1,q = beg;
    for (;b[q]&0x80000000 ; --q)
        a[beg-q]=b[q];
    memcpy(a+beg-q,b,sizeof(unsigned)*(q+1));
    free(BUF);
    return;
}

// int cmp( const void *a , const void *b ) { 
//     return *(float *)a > *(float *)b ? 1 : -1; 
// }
float* sortTwoArray(float* myData, float* sortData, int my_cnt, int sort_cnt, bool isRight) {
    // printf("my_cnt = %d, sort_cnt = %d, rank = %d\n", my_cnt, sort_cnt, rank);
    // float *tmp = (float*) malloc((my_cnt+1) * sizeof(float));
    if(isRight) {
        if(myData[0] >= sortData[sort_cnt-1])
            return myData;
        int my_idx = my_cnt-1;
        int sort_idx = sort_cnt-1;
        int fill_idx = my_cnt-1;
        while(fill_idx >= 0) {
            if(my_idx >= 0 && sort_idx >= 0) {
                if(myData[my_idx] > sortData[sort_idx]) {
                    tmp[fill_idx--] = myData[my_idx--];
                }
                else {
                    tmp[fill_idx--] = sortData[sort_idx--];
                }
            } else if(my_idx >= 0 && sort_idx < 0) {
                tmp[fill_idx--] = myData[my_idx--];
            } else if(my_idx < 0 && sort_idx >= 0) {
                tmp[fill_idx--] = sortData[sort_idx--];
            }
        }
    } else {
        if(myData[my_cnt-1] <= sortData[0])
            return myData;
        int my_idx = 0;
        int sort_idx = 0;
        int fill_idx = 0;
        while(fill_idx < my_cnt) {
            if(my_idx < my_cnt && sort_idx < sort_cnt) {
                if(myData[my_idx] < sortData[sort_idx]) {
                    tmp[fill_idx++] = myData[my_idx++];
                }
                else {
                    tmp[fill_idx++] = sortData[sort_idx++];
                }
            } else if(my_idx < my_cnt && sort_idx >= sort_cnt) {
                tmp[fill_idx++] = myData[my_idx++];
            } else if(my_idx >= my_cnt && sort_idx < sort_cnt) {
                tmp[fill_idx++] = sortData[sort_idx++];
            }
        }
    }
    float* temp = myData;
    myData = tmp;
    tmp = temp;
    return myData;
}

int main(int argc, char *argv[]) {
    const int N = atoi(argv[1]);
    // parameters
    int rank, size;
    int dataPerProc, leftProcDataCount, rightProcDataCount;
    float *proc_data;
    float *buffer;
    MPI_File fin, fout;
    MPI_Comm my_comm = MPI_COMM_WORLD;
    MPI_Group oldGroup, newGroup;

    // MPI init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(size > N) {
        MPI_Comm_group(my_comm, &oldGroup);
        int ranges[][3] = {{ N, size-1, 1 }};
		MPI_Group_range_excl(oldGroup, 1, ranges, &newGroup);
        MPI_Comm_create(my_comm , newGroup , &my_comm);
        size = N;
    }
    if(my_comm == MPI_COMM_NULL){
        MPI_Finalize();
        exit(0);
    }

    // calculate every proccess's data;
    dataPerProc = N/size;
    int leftData = N%size;
    int offset = 0;
    for(int i=0; i<rank; i++) {
        if(i < leftData) {
            offset += (dataPerProc + 1);
        }
        else {
            offset += dataPerProc;
        }
    }
    // neighbor process's data count
    if(rank - 1 < leftData) leftProcDataCount = dataPerProc + 1;
    else leftProcDataCount = dataPerProc;
    if(rank + 1 < leftData) rightProcDataCount = dataPerProc + 1;
    else rightProcDataCount = dataPerProc;

    // how many data in {rank} process
    if(rank < leftData) dataPerProc += 1;
    proc_data = (float*) malloc(dataPerProc * sizeof(float));
    buffer = (float*) malloc((dataPerProc+1) * sizeof(float));
    tmp = (float*) malloc((dataPerProc+1) * sizeof(float));
    // open file and get the data
    MPI_File_open(my_comm, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fin);
    MPI_File_read_at(fin, offset * sizeof(float), proc_data, dataPerProc , MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fin);
    // pre-sort local data
    // qsort(proc_data, dataPerProc, sizeof(float), cmp);
    radix_sort(proc_data, dataPerProc);
    // printf("rank = %d, data num = %d, rightnum = %d, leftnum = %d\n", rank, dataPerProc, rightProcDataCount, leftProcDataCount);
    int sortTimes = 0;
    while(sortTimes <= size) {
        // even phase
        if(rank%2 == 0) { // even index
            if(rank != size-1) {
                MPI_Sendrecv(&proc_data[dataPerProc-1], 1, MPI_FLOAT, rank+1, EvenPhase, buffer, 1, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                // MPI_Send(&proc_data[dataPerProc-1], 1, MPI_FLOAT, rank+1, EvenPhase, my_comm);
                // MPI_Recv(buffer, 1, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                if(proc_data[dataPerProc-1] > buffer[0]) {
                    // printf("rank = %d, proc data = %f, right proc data = %f\n", rank, proc_data[0], buffer[0]);
                    MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank+1, EvenPhase, buffer, rightProcDataCount, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                    // MPI_Send(proc_data, dataPerProc, MPI_FLOAT, rank+1, EvenPhase, my_comm);
                    // MPI_Recv(buffer, rightProcDataCount, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                    proc_data = sortTwoArray(proc_data, buffer, dataPerProc, rightProcDataCount, 0);
                    // printf("rank = %d, proc data = %f\n", rank, proc_data[0]);
                }
            }
        }
        else { // odd index
            MPI_Sendrecv(&proc_data[0], 1, MPI_FLOAT, rank-1, EvenPhase, buffer, 1, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
            // MPI_Recv(buffer, 1, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
            // MPI_Send(&proc_data[0], 1, MPI_FLOAT, rank-1, EvenPhase, my_comm);
            if(proc_data[0] < buffer[0]) {
                // printf("rank = %d, proc data = %f, left proc data = %f\n", rank, proc_data[0], buffer[0]);
                MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank-1, EvenPhase, buffer, leftProcDataCount, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                // MPI_Send(proc_data, dataPerProc, MPI_FLOAT, rank-1, EvenPhase, my_comm);
                // MPI_Recv(buffer, leftProcDataCount, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                proc_data = sortTwoArray(proc_data, buffer, dataPerProc, leftProcDataCount, 1);
                // printf("rank = %d, proc data = %f\n", rank, proc_data[0]);
            }
        }
        // odd phase
        if(rank%2 == 0) { // even index
            if(rank != 0) {
                MPI_Sendrecv(&proc_data[0], 1, MPI_FLOAT, rank-1, EvenPhase, buffer, 1, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                // MPI_Send(&proc_data[0], 1, MPI_FLOAT, rank-1, OddPhase, my_comm);
                // MPI_Recv(buffer, 1, MPI_FLOAT, rank-1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                if(proc_data[0] < buffer[0]) {
                    // printf("rank = %d, proc data = %f, left proc data = %f\n", rank, proc_data[0], buffer[0]);
                    MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank-1, EvenPhase, buffer, leftProcDataCount, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                    // MPI_Send(proc_data, dataPerProc, MPI_FLOAT, rank-1, OddPhase, my_comm);
                    // MPI_Recv(buffer, leftProcDataCount, MPI_FLOAT, rank-1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                    proc_data = sortTwoArray(proc_data, buffer, dataPerProc, leftProcDataCount, 1);
                    // printf("rank = %d, proc data = %f\n", rank, proc_data[0]);
                }
            }
        }
        else { // odd index
            if(rank != size-1) {
                MPI_Sendrecv(&proc_data[dataPerProc-1], 1, MPI_FLOAT, rank+1, EvenPhase, buffer, 1, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                // MPI_Recv(buffer, 1, MPI_FLOAT, rank+1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                // MPI_Send(&proc_data[dataPerProc-1], 1, MPI_FLOAT, rank+1, OddPhase, my_comm);
                if(proc_data[dataPerProc-1] > buffer[0]) {
                    // printf("rank = %d, proc data = %f, right proc data = %f\n", rank, proc_data[0], buffer[0]);
                    MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank+1, EvenPhase, buffer, rightProcDataCount, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                    // MPI_Send(proc_data, dataPerProc, MPI_FLOAT, rank+1, OddPhase, my_comm);
                    // MPI_Recv(buffer, rightProcDataCount, MPI_FLOAT, rank+1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                    proc_data = sortTwoArray(proc_data, buffer, dataPerProc, rightProcDataCount, 0);
                    // printf("rank = %d, proc data = %f\n", rank, proc_data[0]);
                }
            }
        }
        sortTimes += 2;
    }
    // printf("rank = %d, data = %f\n", rank, proc_data[0]);
    MPI_File_open(my_comm, argv[3], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fout);
    MPI_File_write_at(fout, sizeof(float) * offset, proc_data, dataPerProc, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fout);
    MPI_Finalize();
    free(buffer);
    free(proc_data);
    free(tmp);
    return 0;
}
