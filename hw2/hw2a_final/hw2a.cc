#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
#include <emmintrin.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <assert.h>

pthread_mutex_t mutex_use;
pthread_mutex_t mutex_image;

int* image;

typedef struct _myData{
    int iters;
    double left, right;
    double upper, lower;
    int width, height;
    int *now_use, *image;
} myData;

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

void* mandelbrot(void* void_data) {
    myData* data = (myData*) void_data;

    while(1) {
        int rowUse;

        // Get the work batch
        pthread_mutex_lock(&mutex_use);
        rowUse = *(data->now_use);
        *(data->now_use) += 1;
        pthread_mutex_unlock(&mutex_use);

        if(rowUse >= data->height) break;

        bool done[2] = {false, false};
        int repeats[2] = {0,0};
        int rem_point[2] = {0, 0};
        int rowPoint = 0;
        __m128d x, y, x0, y0;
        __m128d length_squared;
        __m128d x_square, y_square;
        __m128d NUM2;
        double range = (data->right - data->left) / data->width;


        // Initialize vec
        x[0] = 0;
        x[1] = 0;
        y0[0] = rowUse * ((data->upper - data->lower) / data->height) + data->lower;
        y0[1] = rowUse * ((data->upper - data->lower) / data->height) + data->lower;
        x0[0] = rowPoint * range + data->left;
        rem_point[0] = rowPoint;
        rowPoint += 1;
        x0[1] = rowPoint * range + data->left;
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

        while(rowPoint <= data->width) {
            y = _mm_add_pd(y0, _mm_mul_pd(NUM2, _mm_mul_pd(x, y)));
            x = _mm_add_pd(_mm_sub_pd(x_square, y_square), x0);
            x_square = _mm_mul_pd(x, x);
            y_square = _mm_mul_pd(y, y);
            length_squared = _mm_add_pd(x_square, y_square);
            repeats[0] += 1;
            repeats[1] += 1;

            if( (length_squared[0] >= 4.0 || repeats[0] >= data->iters) && !done[0]) { 
                // Draw image
                // if(rowUse > 1647) {
                //     printf("RowUse0 = %drem_point0 = %d\n", rowUse, rem_point[1]);
                // }
                image[rowUse * data->width + rem_point[0]] = repeats[0];
                // reset value
                x0[0] = rowPoint * range + data->left;
                rem_point[0] = rowPoint;
                x[0] = 0;
                y[0] = 0;
                length_squared[0] = 0;
                x_square[0] = 0;
                y_square[0] = 0;
                repeats[0] = 0;
                if(rowPoint >= data->width) done[0] = true;
                rowPoint += 1;
            }
            if( (length_squared[1] >= 4.0 || repeats[1] >= data->iters) && !done[1]) { 
                // Draw image
                // if(rowUse > 1647) {
                //     printf("RowUse1 = %drem_point1 = %d\n", rowUse, rem_point[1]);
                // }
                image[rowUse * data->width + rem_point[1]] = repeats[1];
                
                // reset value
                x0[1] = rowPoint * range + data->left;
                rem_point[1] = rowPoint;
                x[1] = 0;
                y[1] = 0;
                length_squared[1] = 0;
                x_square[1] = 0;
                y_square[1] = 0;
                repeats[1] = 0;
                if(rowPoint >= data->width) done[1] = true;
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
            while (no_vec_repeats < data->iters && no_vec_length_squared < 4) {
                no_vec_y = 2 * no_vec_x * no_vec_y + no_vec_y0;
                no_vec_x = no_vec_x_square - no_vec_y_square + no_vec_x0;
                no_vec_x_square = no_vec_x * no_vec_x;
                no_vec_y_square = no_vec_y * no_vec_y;
                no_vec_length_squared = no_vec_x_square + no_vec_y_square;
                ++ no_vec_repeats;
            }
            image[rowUse * data->width + rem_point[0]] = no_vec_repeats;
        }
        if(!done[1]) {
            int no_vec_repeats = repeats[1];
            double no_vec_y0 = y0[1];
            double no_vec_x0 = x0[1];
            double no_vec_x = x[1];
            double no_vec_y = y[1];
            double no_vec_length_squared = length_squared[1];
            double no_vec_x_square = x_square[1];
            double no_vec_y_square = y_square[1];
            while (no_vec_repeats < data->iters && no_vec_length_squared < 4) {
                no_vec_y = 2 * no_vec_x * no_vec_y + no_vec_y0;
                no_vec_x = no_vec_x_square - no_vec_y_square + no_vec_x0;
                no_vec_x_square = no_vec_x * no_vec_x;
                no_vec_y_square = no_vec_y * no_vec_y;
                no_vec_length_squared = no_vec_x_square + no_vec_y_square;
                ++ no_vec_repeats;
            }
            image[rowUse * data->width + rem_point[1]] = no_vec_repeats;
        }
    }
    pthread_exit(NULL);
}

int main(int argc, char** argv) {
    int num_threads;
    /* detect how many CPUs are available */
    cpu_set_t cpu_set;
    sched_getaffinity(0, sizeof(cpu_set), &cpu_set);
    num_threads = CPU_COUNT(&cpu_set);
    printf("%d cpus available\n", CPU_COUNT(&cpu_set));
    pthread_t threads[num_threads];
    int ID[num_threads];

    /* argument parsing */
    assert(argc == 9);
    const char* filename = argv[1];
    int iters = strtol(argv[2], 0, 10);
    double left = strtod(argv[3], 0);
    double right = strtod(argv[4], 0);
    double lower = strtod(argv[5], 0);
    double upper = strtod(argv[6], 0);
    int width = strtol(argv[7], 0, 10);
    int height = strtol(argv[8], 0, 10);

    /* allocate memory for image */
    image = (int*)malloc(width * height * sizeof(int));
    assert(image);
    myData* data = (myData*) malloc (sizeof(myData));
    data->iters = iters;
    data->left = left;
    data->right = right;
    data->lower = lower;
    data->upper = upper;
    data->width = width;
    data->height = height;
    data->now_use = (int*) malloc(sizeof(int));
    *(data->now_use) = 0;
    pthread_mutex_init (&(mutex_use), NULL);
    pthread_mutex_init (&(mutex_image), NULL);
    int t;
    for (t = 0; t < num_threads; t++) {
        ID[t] = t;
        pthread_create(&threads[t], NULL, mandelbrot, (void*)data);
    }
    for (int i=0; i<num_threads; i++) {
		pthread_join(threads[i], NULL);
	}


    /* draw and cleanup */
    write_png(filename, iters, width, height, image);
    pthread_mutex_destroy(&(mutex_use));
    pthread_mutex_destroy(&(mutex_image));
    free(image);
    free(data->now_use);
    free(data);
}
