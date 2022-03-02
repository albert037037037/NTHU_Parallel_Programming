# NTHU_Parallel_Programming
NTHU CS542200 Parallel Programming

## HW1 Odd-Even sort
+  Implement a parallel version of odd-even sort under given restrictions. The goal is to optimize the performance and reduce he total execution time
  + A process can sort or perform any computations on its own local elements.
  + For the odd-even sorting phases, a process can only exchange its local elements with its neighbor processes. Note that the neighbor relationship is not circular
  + However, any communication method, including collective communications (i.e., broadcast, gather, scatter, etc.), are allowed for initialization before the sorting begins, or termination checking.

## HW2 Mandelbrot set
+ Parallized the sequential Mandelbrot Set program by implementing the following two versions:
  1. pthread: Single node shared memory programming using Pthread.
    + This program only needs to be run on a single node.
  2. hybrid: Multi-node hybrid parallelism programming using MPI + OpenMP.
    + This program must be run across multiple nodes.
    + MPI processes are used to balance tasks among nodes, and OpenMP threads are used to perform computations.
    + Pthread library could also be used to create additional threads for communications.
+ Requirements:
  + No mathematical optimization is permitted. That means the computations must be performed on each and every pixel.
 ## HW3 All-Pairs Shortest Path
 + Implement 3 versions of programs that solve the all-pairs shortest path problem.
  + CPU version
    + Use threading to parallelize the computation in your program.
  + Single-GPU version
    + Use Blocked Floyd-Warshall Algorithm.
    + Should be optimized to get the performance points.
  + Multi-GPU version
    + Use Blocked Floyd-Warshall Algorithm.
    + Must use 2 GPUs. Single GPU version is not accepted and will get 0 for correctness and performance score.
 ## HW4 MapReduce
 + Implement a parallel program that mimics the data locality-aware scheduling policy and the functional level programming model of MapReduce
 + Implement the parallel program using MPI and Pthread library.
 + Implement a WordCount sample code to demonstrate your implementation
