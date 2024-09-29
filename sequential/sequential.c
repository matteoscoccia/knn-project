#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

typedef struct {
    float x, y, z;
} Point;

typedef struct {
    float neighbor_distance;
    int neighbor_index;
} NeighborDistance;

// Generate random points in a 3D space
void build_dataset(Point *points, int n) {
    for (int i = 0; i < n; i++) {
        points[i].x = ((float) rand() / RAND_MAX) * 100;
        points[i].y = ((float) rand() / RAND_MAX) * 100;
        points[i].z = ((float) rand() / RAND_MAX) * 100;
    }
}

float euclidean_distance(Point p1, Point p2) {
    return sqrt((p2.x - p1.x) * (p2.x - p1.x) + 
                (p2.y - p1.y) * (p2.y - p1.y) + 
                (p2.z - p1.z) * (p2.z - p1.z));
}


int compare_distances(const void *a, const void *b) {
    NeighborDistance *d1 = (NeighborDistance *)a;
    NeighborDistance *d2 = (NeighborDistance *)b;
    return (d1->neighbor_distance < d2->neighbor_distance) ? -1 : 1;
}

// Find the k-nearest neighbors for each point
void find_k_nearest_neighbors(FILE *f, Point *points, int n, int k) {
    for (int i = 0; i < n; i++) {
        NeighborDistance distances[n - 1];
        int count = 0;

        // Compute all the distances from points[i] to all other points
        for (int j = 0; j < n; j++) {
            if (i != j) { // Skip the distance with itself
                distances[count].neighbor_distance = euclidean_distance(points[i], points[j]);
                distances[count].neighbor_index = j;
                count++;
            }
        }

        // Sort the distances
        qsort(distances, n - 1, sizeof(NeighborDistance), compare_distances);

        //printf("Point %d's %d nearest neighbors:\n", i, k);
        for (int m = 0; m < k; m++) {
            //printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].neighbor_distance);
        }
    }
}

// Computes the time for a knn search iteration
double run_knn(Point *points, int n, int k, FILE *f) {
    clock_t start_time = clock();

    find_k_nearest_neighbors(f, points, n, k);

    clock_t end_time = clock();
    return (double)(end_time - start_time) / CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {

    if (argc < 3 || argc > 4) {
        printf("Please provide a correct number of parameters\n");
        printf("Usage: <number_of_points> <number_of_neighbors> [-ev]\n");
        printf("Use [-ev] if you want to run the algorithm in evaluation mode\n");
        return 1;
    }

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);

    // Validation
    if (n <= 0 || k <= 0 || k >= n) {
        printf("Invalid values: Ensure that n > 0, k > 0, and k < n.\n");
        return 1;
    }

    // Evaluation mode will run the code 5 times
    int evaluation_mode = (argc == 4 && strcmp(argv[3], "-ev") == 0) ? 1 : 0;
    int num_iterations = evaluation_mode ? 5 : 1;

    FILE *f;
    f = fopen("output.csv", "a+");
    char buffer[256]; 
    if (fgets(buffer, sizeof(buffer), f) == NULL) {
        fprintf(f, "Numero di punti (N),Numero di vicini (K),Iterazioni,Tempo di esecuzione(s)\n");
    }
    
    if (f == NULL) {
        printf("Error opening output file\n");
        return 1;
    }

    // Points allocation
    Point *points = (Point *)malloc(n * sizeof(Point));
    if (points == NULL) {
        printf("Memory allocation failure\n");
        return 1;
    }

    double total_time = 0.0;

    for (int iter = 0; iter < num_iterations; iter++) {
        srand(time(NULL) + iter);
        build_dataset(points, n);

        double iteration_time = run_knn(points, n, k, f);
        total_time += iteration_time;

        printf("Iteration %d time: %f seconds\n", iter + 1, iteration_time);
        if (evaluation_mode == 0){
            fprintf(f, "%d,%d,%d,%f\n", n, k, iter + 1, iteration_time);
        }
    }

    // Calculate and print the average time in evaluation mode
    if (evaluation_mode) {
        double average_time = total_time / num_iterations;
        printf("\n%d Points - %d K neighbors - Average time over %d iterations: %f seconds\n", n, k, num_iterations, average_time);
        fprintf(f,"%d,%d,%d,%f\n", n, k, num_iterations, average_time);
    }

    free(points);
    fclose(f);

    return 0;
}
