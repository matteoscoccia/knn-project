#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
    float x, y, z;
} Point;

typedef struct {
    float distance;
    int index;
} Distance;

// Function to generate random points in a 3D space
void generate_points(Point *points, int n) {
    for (int i = 0; i < n; i++) {
        points[i].x = ((float) rand() / RAND_MAX) * 100;
        points[i].y = ((float) rand() / RAND_MAX) * 100;
        points[i].z = ((float) rand() / RAND_MAX) * 100;
    }
}

// Function to calculate Euclidean distance
float euclidean_distance(Point p1, Point p2) {
    return sqrt((p2.x - p1.x) * (p2.x - p1.x) + 
                (p2.y - p1.y) * (p2.y - p1.y) + 
                (p2.z - p1.z) * (p2.z - p1.z));
}

// Comparison function for qsort
int compare_distances(const void *a, const void *b) {
    Distance *d1 = (Distance *)a;
    Distance *d2 = (Distance *)b;
    return (d1->distance < d2->distance) ? -1 : 1;
}

// Function to find k-nearest neighbors
void find_k_nearest_neighbors(FILE *f,Point *points, int n, int k) {
    for (int i = 0; i < n; i++) {
        Distance distances[n - 1];
        int count = 0;

        // Calculate distances from point[i] to all other points
        for (int j = 0; j < n; j++) {
            if (i != j) {
                distances[count].distance = euclidean_distance(points[i], points[j]);
                distances[count].index = j;
                count++;
            }
        }

        // Sort distances to find the k-nearest neighbors
        qsort(distances, n - 1, sizeof(Distance), compare_distances);

        // Output the k-nearest neighbors
        printf("Point %d's %d nearest neighbors:\n", i, k);
        fprintf(f,"Point %d's %d nearest neighbors:\n", i, k);
        for (int m = 0; m < k; m++) {
            printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].index, distances[m].distance);
            fprintf(f,"Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].index, distances[m].distance);
        }
    }
}

int main() {
    FILE *f;
    f = fopen("output.txt", "a+");
    int n = 100;  // Number of points
    int k = 5;    // Number of nearest neighbors to find
    Point *points = (Point *)malloc(n * sizeof(Point));
    
    srand(time(NULL));
    generate_points(points, n);

    clock_t start_time = clock();  // Start measuring time

    // Find the k-nearest neighbors for all points
    find_k_nearest_neighbors(f,points, n, k);

    clock_t end_time = clock();  // End time
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time taken (serial): %f seconds\n", time_taken);
    fprintf(f,"Time taken (serial): %f seconds\n", time_taken);

    free(points);
    fclose(f);
    return 0;
}

