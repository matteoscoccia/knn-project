#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NUM_POINTS 16384
#define MAX_COORDINATE 100.0

typedef struct {
    float x, y, z;
} Point;

void generate_points(Point *points, int n) {
    for (int i = 0; i < n; i++) {
        points[i].x = ((float) rand() / RAND_MAX) * MAX_COORDINATE;
        points[i].y = ((float) rand() / RAND_MAX) * MAX_COORDINATE;
        points[i].z = ((float) rand() / RAND_MAX) * MAX_COORDINATE;
    }
}

void save_points_to_file(const char *filename, Point *points, int n) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file for writing!\n");
        return;
    }

    for (int i = 0; i < n; i++) {
        fprintf(file, "%f %f %f\n", points[i].x, points[i].y, points[i].z);
    }

    fclose(file);
    printf("Successfully saved %d points to %s\n", n, filename);
}

int main() {
    Point points[NUM_POINTS];
    const char *filename = "points_dataset.txt";

    srand(time(NULL));

    // Generate random points
    generate_points(points, NUM_POINTS);

    save_points_to_file(filename, points, NUM_POINTS);

    return 0;
}