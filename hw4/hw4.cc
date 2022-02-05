#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <queue>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unistd.h>
#include <mpi.h>
#include <pthread.h>
#include <chrono>

typedef std::pair<std::string, int> Item;

pthread_mutex_t mutex_num_available;
pthread_mutex_t mutex_task;
pthread_mutex_t mutex_write_file;
pthread_mutex_t mutex_mapper_done;
pthread_cond_t condition_num_available;
pthread_cond_t condition_task;
std::queue<std::pair<int, int>> mapper_tasks; // first: chunkIdx, second: whether has locality
std::vector<std::pair<int, int>> reducer_done_tasks; 
std::queue<std::pair<int, std::pair<int, int>>> mapper_done_task; // first: chunkIdx(taskId), second.first execute time, second.second: how many key-value pair this task generate
std::string input_filename;
int delay;
int num_reducer;
int rank, size;
int chunk_size;

bool cmp (Item a, Item b){
    return a.first < b.first;
}

void* Mapper(void* data) {
    std::chrono::steady_clock::time_point execStart, execEnd;
    int* availableThread = (int*) data;
    std::pair<int, int> task;
    while(true) {            
        pthread_mutex_lock(&mutex_task);
        while(mapper_tasks.empty()) {
            pthread_cond_wait(&condition_task, &mutex_task);
        }
        task = mapper_tasks.front();
        mapper_tasks.pop();
        pthread_mutex_unlock(&mutex_task);

        execStart = std::chrono::steady_clock::now();
        // printf("%d %d\n", task.first, task.second);
        if(task.second == 0) { // recv process chunk has no locality
            sleep(delay);
        }
        std::ifstream input_file(input_filename);
        std::string line;
        std::map<std::string, int> word_count;
        std::string parse;
        std::vector<std::string> words;

        // InputSplit
        for (int i=0; i<(task.first-1)*chunk_size; i++) {
            getline(input_file, line);
        }
        for(int i=0; i<chunk_size; i++) {
            getline(input_file, line);
            size_t pos = 0;
            while ((pos = line.find(" ")) != std::string::npos)
            {
                parse = line.substr(0, pos);
                words.push_back(parse);
                line.erase(0, pos + 1);
            }
            
            if (!line.empty()) {
                words.push_back(line);
            }
        }
        input_file.close();
        // Map // and Partition
        std::vector<std::vector<std::string>> partitionToReducer(num_reducer);
        for (auto word : words) {
            // Partition: 
            // int num = word[0];
            // partitionToReducer[num % num_reducer].push_back(word);

            // Map
            if (word_count.count(word) == 0) {
                word_count[word] = 1;
            }
            else {
                word_count[word]++;
            }
        }

        // Writing Intermediate Files
        // for(int i=0; i<num_reducer; i++) {
        std::string outFilePath = "./intermediate/chunk/";
        std::string chunkId = std::to_string(task.first);
        std::string file = outFilePath + chunkId + ".out";
        std::ofstream writeFile(file);
        for(auto word: word_count) {
            writeFile << word.first << ' ' << word.second << '\n';
        }
        writeFile.close();
        // }

        execEnd = std::chrono::steady_clock::now();
        int execTime = std::chrono::duration_cast<std::chrono::milliseconds>(execEnd-execStart).count();
        pthread_mutex_lock(&mutex_mapper_done);
        mapper_done_task.push(std::make_pair(task.first, std::make_pair(execTime, word_count.size())));
        pthread_mutex_unlock(&mutex_mapper_done);

        pthread_mutex_lock(&mutex_num_available);
        *availableThread += 1;
        pthread_cond_signal(&condition_num_available);
        pthread_mutex_unlock(&mutex_num_available);
    }
    
}



int main(int argc, char **argv)
{
    std::string job_name = std::string(argv[1]);
    num_reducer = std::stoi(argv[2]);
    delay = std::stoi(argv[3]);
    input_filename = std::string(argv[4]);
    chunk_size = std::stoi(argv[5]);
    std::string locality_config_filename = std::string(argv[6]);
    std::string output_dir = std::string(argv[7]);

    int rc;
    rc = MPI_Init(&argc, &argv);
    if (rc != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cpu_set_t cpuset;
    sched_getaffinity(0, sizeof(cpuset), &cpuset);
    int ncpus = CPU_COUNT(&cpuset);

    bool isJobTracker;

    if(rank == size-1) {
        isJobTracker = true;
    } else {
        isJobTracker = false;
    }

    if(isJobTracker) { // JobTracker
        // Read Locality File
        int totalChunk = 0, chunk, locality;
        std::vector<int> chunkNum;
        std::unordered_map<int, int> chunkLocality;
        std::ifstream localityFile(locality_config_filename);
        int total_KV_pairs = 0;
        std::chrono::steady_clock::time_point t1, t2, startProgram, endProgram;
        std::string outputLogFilePath = output_dir + job_name + "-log.out";
        std::ofstream outLogFile(outputLogFilePath);

        outLogFile << time(nullptr) << ",Start_Job" << "," << job_name << "," << size << "," << ncpus << "," << num_reducer << "," << delay << "," << input_filename << "," << chunk_size << "," <<locality_config_filename << "," << output_dir << std::endl;

        startProgram = std::chrono::steady_clock::now();
        while (localityFile >> chunk >> locality) {
            chunkNum.push_back(chunk);
            chunkLocality[chunk] = locality;
            totalChunk += 1;
        }
        localityFile.close();

        // Assign Mappers' Job
        int recv_rank;
        while(!chunkNum.empty()) {
            // std::cout << "JobTracker\n";
            // receive request from other node
            MPI_Recv( &recv_rank , 1 , MPI_INT , MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
            // find the chunk with locality
            int taskChunkIdx = -1;
            int hasLoc = 0;
            for(int i=0; i<chunkNum.size(); i++) {
                if(chunkLocality[chunkNum[i]] == recv_rank) {
                    taskChunkIdx = i;
                    hasLoc = 1;
                    break;
                }
            }
            // if no locality found, get the first one.
            if(taskChunkIdx == -1) taskChunkIdx = 0;
            int send[2];
            send[0] = chunkNum[taskChunkIdx];
            send[1] = hasLoc;
            // printf("chunkNum = %d hasLoc = %d\n", send[0], send[1]);
            outLogFile << time(nullptr) << "," << "Dispatch_MapTask" << "," << send[0] << "," << recv_rank << std::endl;
            MPI_Send(&send, 2, MPI_INT, recv_rank, 1, MPI_COMM_WORLD);
            chunkNum.erase(chunkNum.begin()+taskChunkIdx);
        }
        // Make every node knows there is no chunk.
        for(int i=0; i<size-1; i++) {
            int endMap[2] = {-1, -1};
            MPI_Recv( &recv_rank , 1 , MPI_INT , MPI_ANY_SOURCE , 0 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
            MPI_Send(&endMap, 2, MPI_INT, recv_rank, 1, MPI_COMM_WORLD);
        }

        // After all nodes' mappers end their job, inform the JobTracker
        for(int i=0; i<totalChunk; i++) {
            int info[3];
            MPI_Recv( &info , 3 , MPI_INT , MPI_ANY_SOURCE , 2 , MPI_COMM_WORLD , MPI_STATUS_IGNORE);
            total_KV_pairs += info[2];
            double duration = (double) info[1] / 1000.0;
            outLogFile << time(nullptr) << "," << "Complete_MapTask" << "," << info[0] << "," << duration << std::endl;
            // std::cout << "taskId: " << info[0] << " KV pairs = " << info[1] << std::endl;
        }

        outLogFile << time(nullptr) << "," << "Start_Shuffle" << "," << total_KV_pairs << std::endl;
        t1 = std::chrono::steady_clock::now();
        // Start shuffling
        std::ofstream* intermediateFiles = new std::ofstream[num_reducer];
        for(int i=0; i<num_reducer; i++) {
            std::string path = "./intermediate/" + std::to_string(i+1) + ".out";
            intermediateFiles[i] = std::ofstream(path);
        }

        std::string inputLine;
        for(int i=0; i<totalChunk; i++) {
            std::vector<Item> chunkKV;
            std::string path = "./intermediate/chunk/" + std::to_string(i+1) + ".out";
            std::ifstream chunkfile(path);
            // read chunk file
            while(getline(chunkfile, inputLine)){
                int pos = inputLine.find(" ");
                std::string k = inputLine.substr(0, pos);
                int v = std::stoi(inputLine.substr(pos+1));
                chunkKV.push_back(make_pair(k, v));
            }
            // shuffle
            for(auto kv: chunkKV) {
                int hash = (int)kv.first[0] % num_reducer;
                intermediateFiles[hash] << kv.first << ' ' << kv.second << std::endl;
            }
        }
        t2 = std::chrono::steady_clock::now();
        outLogFile << time(nullptr) << "," << "Finish_Shuffle" << "," << (double)std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000.0 << std::endl;

        for(int i=0; i<num_reducer; i++) {
            intermediateFiles[i].close();
        }

        // Start Reduce phase
        int count = 0;
        while(count < num_reducer) {
            MPI_Recv(&recv_rank, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            count += 1;
            outLogFile << time(nullptr) << "," << "Dispatch_ReduceTask" << "," << count << "," << recv_rank << std::endl;
            MPI_Send(&count, 1, MPI_INT, recv_rank, 4, MPI_COMM_WORLD);
        }

        count = 0;
        while(count < size - 1) {
            count += 1;
            int endReduce = -1;
            MPI_Recv(&recv_rank, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // out << time(nullptr) << "," << "Complete_ReduceTask" << "," << info[0] << "," << info[2] << std::endl;
            MPI_Send(&endReduce, 1, MPI_INT, recv_rank, 4, MPI_COMM_WORLD);
        }
        
        count = num_reducer;
        int msg[2];
        while(count > 0) {
            count -= 1;
            MPI_Recv(&msg, 2, MPI_INT, MPI_ANY_SOURCE, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double duration = (double) msg[1] / 1000.0;
            outLogFile << time(nullptr) << "," << "Complete_ReduceTask" << "," << msg[0] << "," << duration << std::endl;
        }

        endProgram = std::chrono::steady_clock::now();
        outLogFile << time(nullptr) << "," << "FinishJob" << "," << (double)std::chrono::duration_cast<std::chrono::milliseconds>(endProgram-startProgram).count() / 1000.0 << std::endl;
        outLogFile.close();
    }
    else { // TaskTracker
        pthread_mutex_init(&mutex_num_available, NULL);
        pthread_mutex_init(&mutex_write_file, NULL);
        pthread_mutex_init(&mutex_task, NULL);
        pthread_mutex_init(&mutex_mapper_done, NULL);
        pthread_cond_init(&condition_num_available, NULL);
        pthread_cond_init(&condition_task, NULL);
        int* availableThread = new int;
        *availableThread = ncpus - 1; // when no mapper thread using, should be same as mapperThread
        int mapperThread = ncpus - 1;
        pthread_t threads[mapperThread];
        int recvChunkAndLoc[2];

        for(int i=0; i<mapperThread; i++) {
            pthread_create(&threads[i], NULL, &Mapper, (void*)availableThread);
        }

        while (true) {
            pthread_mutex_lock(&mutex_num_available);
            while (*availableThread == 0) {
                pthread_cond_wait(&condition_num_available, &mutex_num_available);
            }
            pthread_mutex_unlock(&mutex_num_available);
            MPI_Send(&rank, 1, MPI_INT, size-1, 0, MPI_COMM_WORLD);
            MPI_Recv(&recvChunkAndLoc, 2, MPI_INT, size-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (recvChunkAndLoc[0] == -1) {
                break;
            }
            pthread_mutex_lock(&mutex_num_available);
            *availableThread -= 1;
            pthread_mutex_unlock(&mutex_num_available);

            std::pair<int, int> m_task(recvChunkAndLoc[0], recvChunkAndLoc[1]);
            pthread_mutex_lock(&mutex_task);
            mapper_tasks.push(m_task);
            pthread_cond_signal(&condition_task);
            pthread_mutex_unlock(&mutex_task);
            // usleep(200);            
        }

        // Send Mapper done message
        while(true) {
            pthread_mutex_lock(&mutex_mapper_done);
            if(*availableThread == mapperThread && mapper_done_task.size() == 0) {
                pthread_mutex_unlock(&mutex_mapper_done);
                break;
            }
            else if(mapper_done_task.size() != 0) {
                std::pair<int, std::pair<int, int>> done = mapper_done_task.front();
                mapper_done_task.pop();
                pthread_mutex_unlock(&mutex_mapper_done);
                int mapper_done_signal[3];
                mapper_done_signal[0] = done.first;
                mapper_done_signal[1] = (done.second).first;
                mapper_done_signal[2] = (done.second).second;
                MPI_Send(&mapper_done_signal, 3, MPI_INT, size-1, 2, MPI_COMM_WORLD);
            }
            else pthread_mutex_unlock(&mutex_mapper_done);
        }

        // Start Reduce phase
        std::chrono::steady_clock::time_point reduceExecStart, reduceExecEnd;
        reduceExecStart = std::chrono::steady_clock::now();
        while(true) {
            int getReduceTaskId;
            MPI_Send(&rank, 1, MPI_INT, size-1, 3, MPI_COMM_WORLD);
            MPI_Recv(&getReduceTaskId, 1, MPI_INT, size-1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if(getReduceTaskId == -1) {
                break;
            }
            // std::cout << getReduceTaskId << std::endl;
            std::string filePath = "./intermediate/" + std::to_string(getReduceTaskId) + ".out";
            std::ifstream interFile(filePath);
            std::string input;
            std::vector<Item> KVPair;
            while(getline(interFile, input)) {
                int pos = input.find(" ");
                std::string word = input.substr(0, pos);
                std::string freq_str = input.substr(pos+1);
                if(freq_str == "\n" || freq_str == "") continue;
                int word_freq = std::stoi(freq_str);
                KVPair.push_back(std::make_pair(word, word_freq));
            }

            // Sort the intermediate file
            std::sort(KVPair.begin(), KVPair.end(), cmp);

            // Group
            std::map<std::string, std::vector<int>> groupResult;
            for(auto pair: KVPair) {
                groupResult[pair.first].push_back(pair.second);
            }

            // Reduce
            std::map<std::string, int> reduceResult;
            for(auto group: groupResult) {
                reduceResult[group.first] = 0;
                for(auto v: group.second) {
                    reduceResult[group.first] += v;
                }
            }

            // Output
            std::string outputPath = output_dir + job_name + "-" + std::to_string(getReduceTaskId) + ".out";
            std::ofstream outputFile(outputPath);
            for(auto iter: reduceResult) {
                outputFile << iter.first << " " << iter.second << std::endl;
            }
            reduceExecEnd = std::chrono::steady_clock::now();
            double reduceExecTime = std::chrono::duration_cast<std::chrono::milliseconds>(reduceExecEnd-reduceExecStart).count();
            reducer_done_tasks.push_back(std::make_pair(getReduceTaskId, reduceExecTime));
        }

        while(!reducer_done_tasks.empty()) {
            std::pair <int, double> doneTask = reducer_done_tasks.back();
            reducer_done_tasks.pop_back();
            int sendLog[2];
            sendLog[0] = doneTask.first;
            sendLog[1] = doneTask.second;
            MPI_Send(sendLog, 2, MPI_INT, size-1, 5, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}
