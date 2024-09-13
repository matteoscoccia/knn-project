#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Define the structure for a 3D point
typedef struct {
    float x, y, z;
} Point;

// Define the structure for storing the distance and the index of a point
typedef struct {
    float neighbor_distance;
    int neighbor_index;
} NeighborDistance;

// Function to generate random points in a 3D space
void generate_points(Point *points, int n) {
    for (int i = 0; i < n; i++) {
        points[i].x = ((float) rand() / RAND_MAX) * 100;
        points[i].y = ((float) rand() / RAND_MAX) * 100;
        points[i].z = ((float) rand() / RAND_MAX) * 100;
    }
}

// Function to calculate Euclidean distance between two 3D points
float euclidean_distance(Point p1, Point p2) {
    return sqrt((p2.x - p1.x) * (p2.x - p1.x) + 
                (p2.y - p1.y) * (p2.y - p1.y) + 
                (p2.z - p1.z) * (p2.z - p1.z));
}

// Comparison function for sorting distances using qsort
int compare_distances(const void *a, const void *b) {
    NeighborDistance *d1 = (NeighborDistance *)a;
    NeighborDistance *d2 = (NeighborDistance *)b;
    return (d1->neighbor_distance < d2->neighbor_distance) ? -1 : 1;
}

void insertion_sort(NeighborDistance *distances, int size) {
    for (int i = 1; i < size; i++) {
        NeighborDistance key = distances[i];
        int j = i - 1;

        // Shift elements of distances[0..i-1], that are greater than key, to one position ahead of their current position
        while (j >= 0 && distances[j].neighbor_distance > key.neighbor_distance) {
            distances[j + 1] = distances[j];
            j = j - 1;
        }
        distances[j + 1] = key;
    }
}

// Function to find the k-nearest neighbors for each point
void find_k_nearest_neighbors(FILE *f, Point *points, int n, int k) {
    for (int i = 0; i < n; i++) {
        NeighborDistance distances[n - 1];
        int count = 0;

        // Calculate distances from points[i] to all other points
        for (int j = 0; j < n; j++) {
            if (i != j) {
                distances[count].neighbor_distance = euclidean_distance(points[i], points[j]);
                distances[count].neighbor_index = j;
                count++;
            }
        }

        // Sort the distances to find the k-nearest neighbors
        qsort(distances, n - 1, sizeof(NeighborDistance), compare_distances);
        //insertion_sort(distances, n - 1);

        // Output the k-nearest neighbors for the current point
        printf("Point %d's %d nearest neighbors:\n", i, k);
        fprintf(f, "Point %d's %d nearest neighbors:\n", i, k);
        for (int m = 0; m < k; m++) {
            printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].neighbor_distance);
            fprintf(f, "Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].neighbor_distance);
        }
    }
}

int main(int argc, char *argv[]) {
    // Check if the number of command-line arguments is correct
    if (argc != 3) {
        printf("Please provide the correct parameters. Usage: %s <number_of_points> <number_of_neighbors>\n", argv[0]);
        return 1;
    }

    // Parse the number of points (n) and number of neighbors (k) from command-line arguments
    int n = atoi(argv[1]);
    int k = atoi(argv[2]);

    // Validate the input values
    if (n <= 0 || k <= 0 || k >= n) {
        printf("Invalid values: Ensure that n > 0, k > 0, and k < n.\n");
        return 1;
    }

    FILE *f;
    f = fopen("output.txt", "a+");
    if (f == NULL) {
        printf("Error opening output file!\n");
        return 1;
    }

    // Allocate memory for the points
    Point *points = (Point *)malloc(n * sizeof(Point));
    if (points == NULL) {
        printf("Memory allocation failed!\n");
        return 1;
    }

    srand(time(NULL));  // Seed for random number generation
    generate_points(points, n);  // Generate random points

    clock_t start_time = clock();  // Start measuring time

    // Find the k-nearest neighbors for all points
    find_k_nearest_neighbors(f, points, n, k);

    clock_t end_time = clock();  // End time
    double time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time taken (serial): %f seconds\n", time_taken);
    fprintf(f, "Time taken (serial): %f seconds\n", time_taken);

    // Free the allocated memory and close the file
    free(points);
    fclose(f);

    return 0;
}
