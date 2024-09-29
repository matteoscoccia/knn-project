#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

typedef struct {
    float x;
    float y;
    float z;
} Point;

typedef struct {
    float neighbor_distance;
    int neighbor_index;
} Distance;

// Function to generate random points in a 3D space
void generate_points(Point *points, int num_points) {
    for (int i = 0; i < num_points; i++) {
        points[i].x = ((float) rand() / RAND_MAX) * 100;
        points[i].y = ((float) rand() / RAND_MAX) * 100;
        points[i].z = ((float) rand() / RAND_MAX) * 100;
    }
}

// Function to calculate Euclidean distance between two points
float euclidean_distance(Point *point1, Point *point2) {
    float x = point1->x - point2->x;
    float y = point1->y - point2->y;
    float z = point1->z - point2->z;
    return sqrt(x * x + y * y + z * z);
}

// Insert a neighbor into the K nearest neighbors list if it is closer
void insert_distance_if_needed(Distance *nearest_neighbors, int K, Distance new_neighbor) {
    if (new_neighbor.neighbor_distance < nearest_neighbors[K-1].neighbor_distance) {
        nearest_neighbors[K-1] = new_neighbor;

        // Move the new element to its correct position to maintain sorted order
        for (int i = K - 2; i >= 0; i--) {
            if (nearest_neighbors[i].neighbor_distance > nearest_neighbors[i + 1].neighbor_distance) {
                Distance temp = nearest_neighbors[i];
                nearest_neighbors[i] = nearest_neighbors[i + 1];
                nearest_neighbors[i + 1] = temp;
            } else {
                break;
            }
        }
    }
}

// Calculate distances and retain K nearest neighbors
void calculate_and_retain_k_nearest(Point *points, Distance *results, int N, int K) {
    for (int i = 0; i < N; i++) {
        Distance *nearest_neighbors = &results[i * K];

        for (int k = 0; k < K; k++) {
            nearest_neighbors[k].neighbor_distance = INFINITY;
            nearest_neighbors[k].neighbor_index = -1;
        }

        // Calculate the distance of the current point to all other points
        for (int j = 0; j < N; j++) {
            if (i != j) {
                Distance new_neighbor;
                new_neighbor.neighbor_distance = euclidean_distance(&points[i], &points[j]);
                new_neighbor.neighbor_index = j;

                insert_distance_if_needed(nearest_neighbors, K, new_neighbor);
            }
        }
    }
}

// Run KNN and measure execution time
double run_knn(Point *points, Distance *results, int N, int K) {
    clock_t start_time = clock();

    calculate_and_retain_k_nearest(points, results, N, K);

    clock_t end_time = clock();
    return (double)(end_time - start_time) / CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {
    if (argc < 3 || argc > 4) {
        printf("Please provide the correct number of parameters\n");
        printf("Usage: <number_of_points> <number_of_neighbors> [-ev]\n");
        return 1;
    }

    int N = atoi(argv[1]);
    int K = atoi(argv[2]);

    // Check if K is valid
    if (K >= N) {
        printf("K must be less than N.\n");
        return 1;
    }

    // Check if evaluation mode is enabled (with -ev parameter)
    int evaluation_mode = (argc == 4 && strcmp(argv[3], "-ev") == 0) ? 1 : 0;
    int num_iterations = evaluation_mode ? 5 : 1;

    FILE *f;
    f = fopen("output_sequential_improved.csv", "a+");
    char buffer[256]; 
    if (fgets(buffer, sizeof(buffer), f) == NULL) {
        fprintf(f, "Numero di punti (N),Numero di vicini (K),Iterazioni,Tempo di esecuzione(s)\n");
    }
    
    if (f == NULL) {
        printf("Error opening output file\n");
        return 1;
    }


    Point *points = (Point *)malloc(N * sizeof(Point));
    Distance *results = (Distance *)malloc(N * K * sizeof(Distance));

    double total_time = 0.0;

    if(evaluation_mode){
        printf("Running in evaluation mode\n");
    }

    for (int iter = 0; iter < num_iterations; iter++) {
        srand(time(NULL) + iter);
        generate_points(points, N);  // Generate random points

        double iteration_time = run_knn(points, results, N, K);
        total_time += iteration_time;

        printf("Iteration %d time: %f seconds\n", iter + 1, iteration_time);

        if (!evaluation_mode) {
            /*printf("\nResults for iteration %d:\n", iter + 1);
            for (int i = 0; i < N; i++) {
                printf("Point %d's %d nearest neighbors:\n", i, K);
                for (int m = 0; m < K; m++) {
                    printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, results[i * K + m].neighbor_index, results[i * K + m].neighbor_distance);
                }
            }*/
            fprintf(f, "%d,%d,%d,%f\n", N, K, iter + 1, iteration_time);
        }
    }

    
    if (evaluation_mode) {
        double average_time = total_time / num_iterations;
        printf("\n%d Points - %d K neighbors - Average time over %d iterations: %f seconds\n", N, K, num_iterations, average_time);
        fprintf(f,"%d,%d,%d,%f\n", N, K, num_iterations, average_time);
    }

    // Free allocated memory
    free(points);
    free(results);

    return 0;
}
