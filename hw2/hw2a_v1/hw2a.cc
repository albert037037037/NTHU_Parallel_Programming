#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define PNG_NO_SETJMP
#include <sched.h>
#include <assert.h>
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

void* cal_mandelbrot(void* void_data) {
    myData* data = (myData*) void_data;

    while(*(data->now_use) < data->width * data->height) {
        pthread_mutex_lock(&mutex_use);
        int startj = *(data->now_use) % data->width;
        int starti = *(data->now_use) / data->width;
        if(*(data->now_use) + 50 >= data->width * data->height){
            *(data->now_use) = data->width * data->height;
        }
        else {
            *(data->now_use) += 50;
        }
        int endj = *(data->now_use) % data->width;
        int endi = *(data->now_use) / data->width;
        pthread_mutex_unlock(&mutex_use);

        double y0, x0;
        if(starti == endi) {
            y0 = starti * ((data->upper - data->lower) / data->height) + data->lower;
            for(int j=startj; j<endj; j++) {
                x0 = j * ((data->right - data->left) / data->width) + data->left;
                int repeats = 0;
                double x = 0;
                double y = 0;
                double length_squared = 0;
                while (repeats < data->iters && length_squared < 4) {
                    double temp = x * x - y * y + x0;
                    y = 2 * x * y + y0;
                    x = temp;
                    length_squared = x * x + y * y;
                    ++repeats;
                }
                image[starti * data->width + j] = repeats;
            }
        }
        else {
            y0 = starti * ((data->upper - data->lower) / data->height) + data->lower;
            for(int j=startj; j<data->width; j++) {
                x0 = j * ((data->right - data->left) / data->width) + data->left;
                int repeats = 0;
                double x = 0;
                double y = 0;
                double length_squared = 0;
                while (repeats < data->iters && length_squared < 4) {
                    double temp = x * x - y * y + x0;
                    y = 2 * x * y + y0;
                    x = temp;
                    length_squared = x * x + y * y;
                    ++repeats;
                }
                image[starti * data->width + j] = repeats;
            }
            y0 = endi * ((data->upper - data->lower) / data->height) + data->lower;
            for(int j=0; j<endj; j++) {
                x0 = j * ((data->right - data->left) / data->width) + data->left;
                int repeats = 0;
                double x = 0;
                double y = 0;
                double length_squared = 0;
                while (repeats < data->iters && length_squared < 4) {
                    double temp = x * x - y * y + x0;
                    y = 2 * x * y + y0;
                    x = temp;
                    length_squared = x * x + y * y;
                    ++repeats;
                }
                image[endi * data->width + j] = repeats;
            }
        }
    }
    return NULL;
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
        pthread_create(&threads[t], NULL, cal_mandelbrot, (void*)data);
    }
    for (int i=0; i<num_threads; i++) {
		pthread_join(threads[i], NULL);
	}
    // for(int i=0; i<height; i++) {
    //     for(int j=0; j<width; j++) {
    //         printf("%d ", image[i*width+j]);
    //     }
    //     printf("\n");
    // }
    /* draw and cleanup */
    write_png(filename, iters, width, height, image);
    pthread_mutex_destroy(&(mutex_use));
    pthread_mutex_destroy(&(mutex_image));
    free(image);
}
