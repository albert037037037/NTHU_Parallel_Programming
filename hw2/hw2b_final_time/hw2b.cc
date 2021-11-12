#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <emmintrin.h>
#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include <chrono>

void write_png(const char* filename, int iters, int width, int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_ptr, 0, PNG_NO_FILTERS);
    png_write_info(png_ptr, info_ptr);
    png_set_compression_level(png_ptr, 1);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            png_bytep color = row + x * 3;
            if (p != iters) {
                if (p & 16) {
                    color[0] = 240;
                    color[1] = color[2] = p % 16 * 16;
                } else {
                    color[0] = p % 16 * 16;
                }
            }
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

int main(int argc, char** argv) {
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    printf("%d cpus available\n", CPU_COUNT(&cpu_set));
    int num_thread = CPU_COUNT(&cpu_set);

    // Initialize MPI settings
    MPI_Init(&argc, &argv);
    int rank, size, omp_threads, omp_thread;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::chrono::steady_clock::time_point t1, t2;

    /* argument parsing */
    assert(argc == 9);
    t1 = std::chrono::steady_clock::now();
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);

    

    /* allocate memory for image */
    int* image = (int*)malloc(width * height * sizeof(int));
    memset(image, 0, width * height * sizeof(int));
    assert(image);

#pragma omp parallel for schedule(dynamic) num_threads(num_thread) shared(iters, left, right, lower, upper, width, height)
        for(int rowUse=rank; rowUse<height; rowUse+=size) {

            bool done[2] = {false, false};
            int repeats[2] = {0, 0};
            int rem_point[2] = {0, 0};
            int rowPoint = 0;
            __m128d x, y, x0, y0;
            __m128d length_squared;
            __m128d x_square, y_square;
            __m128d NUM2;

            // Initialize vec
            x[0] = 0;
            x[1] = 0;
            y0[0] = rowUse * ((upper - lower) / height) + lower;
            y0[1] = rowUse * ((upper - lower) / height) + lower;
            x0[0] = rowPoint * ((right - left) / width) + left;
            rem_point[0] = rowPoint;
            rowPoint += 1;
            x0[1] = rowPoint * ((right - left) / width) + left;
            rem_point[1] = rowPoint;
            rowPoint += 1;
            y[0] = 0;
            y[1] = 0;
            length_squared[0] = 0;
            length_squared[1] = 0;
            x_square[0] = 0;
            x_square[1] = 0;
            y_square[0] = 0;
            y_square[1] = 0;
            NUM2[0] = 2.0;
            NUM2[1] = 2.0;

            while(rowPoint <= width) {
                y = _mm_add_pd(y0, _mm_mul_pd(NUM2, _mm_mul_pd(x, y)));
                x = _mm_add_pd(_mm_sub_pd(x_square, y_square), x0);
                x_square = _mm_mul_pd(x, x);
                y_square = _mm_mul_pd(y, y);
                length_squared = _mm_add_pd(x_square, y_square);
                repeats[0] += 1;
                repeats[1] += 1;

                if( (length_squared[0] >= 4.0 || repeats[0] >= iters) && !done[0]) { 
                    // Draw image
                    image[rowUse * width + rem_point[0]] = repeats[0];
                    // reset value
                    x0[0] = rowPoint * ((right - left) / width) + left;
                    rem_point[0] = rowPoint;
                    x[0] = 0;
                    y[0] = 0;
                    length_squared[0] = 0;
                    x_square[0] = 0;
                    y_square[0] = 0;
                    repeats[0] = 0;
                    if(rowPoint >= width) done[0] = true;
                    rowPoint += 1;
                }
                if( (length_squared[1] >= 4.0 || repeats[1] >= iters) && !done[1]) { 
                    // Draw image
                    image[rowUse * width + rem_point[1]] = repeats[1];
                    
                    // reset value
                    x0[1] = rowPoint * ((right - left) / width) + left;
                    rem_point[1] = rowPoint;
                    x[1] = 0;
                    y[1] = 0;
                    length_squared[1] = 0;
                    x_square[1] = 0;
                    y_square[1] = 0;
                    repeats[1] = 0;
                    if(rowPoint >= width) done[1] = true;
                    rowPoint += 1;
                }
            }

            if(!done[0]) {
                int no_vec_repeats = repeats[0];
                double no_vec_y0 = y0[0];
                double no_vec_x0 = x0[0];
                double no_vec_x = x[0];
                double no_vec_y = y[0];
                double no_vec_length_squared = length_squared[0];
                double no_vec_x_square = x_square[0];
                double no_vec_y_square = y_square[0];
                while (no_vec_repeats < iters && no_vec_length_squared < 4) {
                    no_vec_y = 2 * no_vec_x * no_vec_y + no_vec_y0;
                    no_vec_x = no_vec_x_square - no_vec_y_square + no_vec_x0;
                    no_vec_x_square = no_vec_x * no_vec_x;
                    no_vec_y_square = no_vec_y * no_vec_y;
                    no_vec_length_squared = no_vec_x_square + no_vec_y_square;
                    ++ no_vec_repeats;
                }
                image[rowUse * width + rem_point[0]] = no_vec_repeats;
            }
            else if(!done[1]) {
                int no_vec_repeats = repeats[1];
                double no_vec_y0 = y0[1];
                double no_vec_x0 = x0[1];
                double no_vec_x = x[1];
                double no_vec_y = y[1];
                double no_vec_length_squared = length_squared[1];
                double no_vec_x_square = x_square[1];
                double no_vec_y_square = y_square[1];
                while (no_vec_repeats < iters && no_vec_length_squared < 4) {
                    no_vec_y = 2 * no_vec_x * no_vec_y + no_vec_y0;
                    no_vec_x = no_vec_x_square - no_vec_y_square + no_vec_x0;
                    no_vec_x_square = no_vec_x * no_vec_x;
                    no_vec_y_square = no_vec_y * no_vec_y;
                    no_vec_length_squared = no_vec_x_square + no_vec_y_square;
                    ++ no_vec_repeats;
                }
                image[rowUse * width + rem_point[1]] = no_vec_repeats;
            }
        }
    


    t2 = std::chrono::steady_clock::now();
    long int compute_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    long int sumCommunicate = 0;
    MPI_Reduce(&compute_time, &sumCommunicate, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0) {
        double avgCommunicate = double(sumCommunicate) / (double)size;
        printf("Compute %.2lf ms.\n", avgCommunicate);
    }
    t1 = std::chrono::steady_clock::now();
    if(rank == 0) MPI_Reduce(MPI_IN_PLACE, image, width*height, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);  
    else MPI_Reduce(image, NULL, width*height, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    t2 = std::chrono::steady_clock::now();
    printf("Communicate %ld ms.\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
    /* draw and cleanup */
    if(rank == 0) {
        t1 = std::chrono::steady_clock::now();
        write_png(filename, iters, width, height, image);
        t2 = std::chrono::steady_clock::now();
        printf("IO %ld ms.\n", std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());
    }
    free(image);
}
