#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

typedef struct {
    float x, y, z;
} Point;

// Contains the distance from the point and the index of the point
typedef struct {
    float neighbor_distance;
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
    return (d1->neighbor_distance < d2->neighbor_distance) ? -1 : 1;
}

// Function to find k-nearest neighbors
void find_k_nearest_neighbors(Distance *local_results, Point *points, int n, int k, int start, int end) {
    for (int i = start; i < end; i++) {
        Distance distances[n - 1];
        int count = 0;

        for (int j = 0; j < n; j++) {
            if (i != j) {
                distances[count].neighbor_distance = euclidean_distance(points[i], points[j]);
                distances[count].neighbor_index = j;
                count++;
            }
        }

        // Sort distances
        qsort(distances, n - 1, sizeof(Distance), compare_distances);

        printf("Point %d's %d nearest neighbors:\n", i, k);
        // Copy the k-nearest neighbors to the local_results array
        for (int m = 0; m < k; m++) {
            local_results[(i - start) * k + m] = distances[m];            
            printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].neighbor_distance);
        }

    }
}

int main(int argc, char** argv) {
    FILE *f;
    f = fopen("output_parallel.txt", "a+");
    if (f == NULL) {
        printf("Error opening output file!\n");
        return 1;
    }

    int rank, size, n = 16384, k = 20;
    Point *points = NULL;
    Distance *local_results = NULL;
    Distance *global_results = NULL;
    double start_time, end_time, local_time, max_time;

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
    if (rank == size - 1) end = n;  // Last process handles remainings

    // Allocate local results array
    local_results = (Distance*)malloc(points_per_process * k * sizeof(Distance));

    start_time = MPI_Wtime();  // Start measuring parallel time

    //TODO dichiarare MPI TYPE DISTANCE

    // Each process finds k-nearest neighbors for its share of points
    find_k_nearest_neighbors(local_results, points, n, k, start, end);


    // Gather results from all processes
    if (rank == 0) {
        global_results = (Distance*)malloc(n * k * sizeof(Distance));
    }
    MPI_Gather(local_results, points_per_process * k * sizeof(Distance), MPI_BYTE,
               global_results, points_per_process * k * sizeof(Distance), MPI_BYTE,
               0, MPI_COMM_WORLD);

    /*if (rank == size - 1) {
        int remaining_points = n - (size - 1) * points_per_process;
        MPI_Gatherv(local_results, remaining_points * k * sizeof(Distance), MPI_BYTE,
                    global_results, NULL, NULL, MPI_BYTE,
                    0, MPI_COMM_WORLD);
    }*/

    end_time = MPI_Wtime();  // End time
    local_time = end_time - start_time;

    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            printf("Point %d's %d nearest neighbors:\n", i, k);
            for (int m = 0; m < k; m++) {
                printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, global_results[i * k + m].neighbor_index, global_results[i * k + m].neighbor_distance);
            }
        }
        printf("Parallel time taken: %f seconds\n", max_time);
        fprintf(f, "Parallel time taken: %f seconds\n", max_time);
        
        printf("\n%d Points - %d K neighbors - Average time over 1 iterations: %f seconds\n", n, k, max_time);
        fprintf(f, "%d Points - %d K neighbors - Average time over 1 iterations: %f seconds\n", n, k, max_time);
        free(global_results);
    }

    free(local_results);
    free(points);
    MPI_Finalize();
    fclose(f);
    return 0;
}
