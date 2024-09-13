#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

typedef struct {
    float x, y, z;
} Point;

//Contiene la distanza dal punto e l'indice del punto
typedef struct {
    float distance;
    int neighbor_index;
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
void find_k_nearest_neighbors( FILE *f,Point *points, int n, int k, int start, int end) {
    for (int i = start; i < end; i++) {
        Distance distances[n - 1];
        int count = 0;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                distances[count].distance = euclidean_distance(points[i], points[j]);
                distances[count].neighbor_index = j;
                count++;
            }
        }

        // Sort distances
        qsort(distances, n - 1, sizeof(Distance), compare_distances);

        // Print the k-nearest neighbors
        printf("Process %d, Point %d's %d nearest neighbors:\n", start, i, k);
        for (int m = 0; m < k; m++) {
            printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].distance);
            fprintf(f,"Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].distance);
        }
    }
}

int main(int argc, char** argv) {
    FILE *f;
    f = fopen("output.txt", "a+");
    int rank, size, n = 100, k = 5;
    Point *points = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Root process generates points
    if (rank == 0) {
        points = (Point*)malloc(n * sizeof(Point));
        srand(time(NULL));
        generate_points(points, n);
    }

    // Broadcast points to all processes
    if (rank != 0) {
        points = (Point*)malloc(n * sizeof(Point));
    }
    MPI_Bcast(points, n * sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Divide the work among processes
    int points_per_process = n / size;
    int start = rank * points_per_process;
    int end = (rank + 1) * points_per_process;
    if (rank == size - 1) end = n;  // Last process handles remainder

    double start_time = MPI_Wtime();  // Start measuring parallel time

    // Each process finds k-nearest neighbors for its share of points
    find_k_nearest_neighbors(f,points, n, k, start, end);

    double end_time = MPI_Wtime();  // End time
    if (rank == 0) {
        printf("Parallel time taken: %f seconds\n", end_time - start_time);
        fprintf(f,"Parallel time taken: %f seconds\n", end_time - start_time);
    }

    free(points);
    MPI_Finalize();
    fclose(f);
    return 0;
}

