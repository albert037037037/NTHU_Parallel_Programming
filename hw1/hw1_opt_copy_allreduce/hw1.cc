#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <boost/sort/spreadsort/float_sort.hpp>

#define EvenPhase 998
#define OddPhase 999
float *tmp;

float* sortTwoArray(float* myData, float* sortData, int my_cnt, int sort_cnt, bool isRight) {
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
    bool modified, all_modified; // to check whether sorting is done
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
    boost::sort::spreadsort::detail::float_sort(proc_data, &proc_data[dataPerProc]);
    // printf("rank = %d, data num = %d, rightnum = %d, leftnum = %d\n", rank, dataPerProc, rightProcDataCount, leftProcDataCount);
    while(true) {
        modified = false;
        // even phase
        if(rank%2 == 0) { // even index
            if(rank != size-1) {
                MPI_Sendrecv(&proc_data[dataPerProc-1], 1, MPI_FLOAT, rank+1, EvenPhase, buffer, 1, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                if(proc_data[dataPerProc-1] > buffer[0]) {
                    modified = true;
                    MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank+1, EvenPhase, buffer, rightProcDataCount, MPI_FLOAT, rank+1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                    proc_data = sortTwoArray(proc_data, buffer, dataPerProc, rightProcDataCount, 0);
                }
            }
        }
        else { // odd index
            MPI_Sendrecv(&proc_data[0], 1, MPI_FLOAT, rank-1, EvenPhase, buffer, 1, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
            if(proc_data[0] < buffer[0]) {
                modified = true;
                MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank-1, EvenPhase, buffer, leftProcDataCount, MPI_FLOAT, rank-1, EvenPhase, my_comm, MPI_STATUS_IGNORE);
                proc_data = sortTwoArray(proc_data, buffer, dataPerProc, leftProcDataCount, 1);
            }
        }
        // odd phase
        if(rank%2 == 0) { // even index
            if(rank != 0) {
                MPI_Sendrecv(&proc_data[0], 1, MPI_FLOAT, rank-1, OddPhase, buffer, 1, MPI_FLOAT, rank-1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                if(proc_data[0] < buffer[0]) {
                    modified = true;
                    MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank-1, OddPhase, buffer, leftProcDataCount, MPI_FLOAT, rank-1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                    proc_data = sortTwoArray(proc_data, buffer, dataPerProc, leftProcDataCount, 1);
                }
            }
        }
        else { // odd index
            if(rank != size-1) {
                MPI_Sendrecv(&proc_data[dataPerProc-1], 1, MPI_FLOAT, rank+1, OddPhase, buffer, 1, MPI_FLOAT, rank+1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                if(proc_data[dataPerProc-1] > buffer[0]) {
                    modified = true;
                    MPI_Sendrecv(proc_data, dataPerProc, MPI_FLOAT, rank+1, OddPhase, buffer, rightProcDataCount, MPI_FLOAT, rank+1, OddPhase, my_comm, MPI_STATUS_IGNORE);
                    proc_data = sortTwoArray(proc_data, buffer, dataPerProc, rightProcDataCount, 0);
                }
            }
        }
        all_modified = false;
        MPI_Allreduce( &modified , &all_modified , 1 , MPI_C_BOOL , MPI_LOR , my_comm);
        if(all_modified == false) break;
    }
    MPI_File_open(my_comm, argv[3], MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fout);
    MPI_File_write_at(fout, sizeof(float) * offset, proc_data, dataPerProc, MPI_FLOAT, MPI_STATUS_IGNORE);
    MPI_File_close(&fout);
    MPI_Finalize();
    free(buffer);
    free(proc_data);
    free(tmp);
    return 0;
}
