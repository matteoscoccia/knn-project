#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <string.h>
#include <time.h>

typedef struct {
    float x, y, z;
} Point;

// Contains the distance from the point and the index of the point
typedef struct {
    float neighbor_distance;
    int neighbor_index;
} Distance;

// Generate random points in a 3D space
void generate_points(Point *points, int n) {
    printf("\n----Generating %d points----\n", n);
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

// Comparison function for qsort
int compare_distances(const void *a, const void *b) {
    Distance *d1 = (Distance *)a;
    Distance *d2 = (Distance *)b;
    return (d1->neighbor_distance < d2->neighbor_distance) ? -1 : 1;
}

// Create MPI datatype for Point structure
void create_mpi_point_type(MPI_Datatype *mpi_point_type) {
    int block_lengths[3] = {1, 1, 1};  // x, y, z are all floats
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Aint offsets[3];

    // Offsets of structure members
    offsets[0] = offsetof(Point, x);
    offsets[1] = offsetof(Point, y);
    offsets[2] = offsetof(Point, z);

    // Create MPI datatype for Point structure
    MPI_Type_create_struct(3, block_lengths, offsets, types, mpi_point_type);
    MPI_Type_commit(mpi_point_type);
}

// Create MPI datatype for Distance structure
void create_mpi_distance_type(MPI_Datatype *mpi_distance_type) {
    int block_lengths[2] = {1, 1};  // 1 float and 1 int
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_INT};
    MPI_Aint offsets[2];

    // Offsets of structure members
    offsets[0] = offsetof(Distance, neighbor_distance);
    offsets[1] = offsetof(Distance, neighbor_index);

    // Create MPI datatype for Distance structure
    MPI_Type_create_struct(2, block_lengths, offsets, types, mpi_distance_type);
    MPI_Type_commit(mpi_distance_type);
}

// Find k-nearest neighbors
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

        //printf("Point %d's %d nearest neighbors:\n", i, k);
        // Copy the k-nearest neighbors to the local_results array
        for (int m = 0; m < k; m++) {
            local_results[(i - start) * k + m] = distances[m];            
            //printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, distances[m].neighbor_index, distances[m].neighbor_distance);
        }
    }
}

int main(int argc, char** argv) {

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

    printf("\nStarting");


    int rank, size;
    Point *points = NULL;
    Distance *local_results = NULL;
    Distance *global_results = NULL;
    double start_time, end_time, local_time, max_time;
    MPI_Datatype mpi_distance_type, mpi_point_type;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("\nSize = %d", size);
    printf("\nN = %d", n);
    printf("\nP = %d", size);
    printf("\nK = %d", k);
    
    FILE *f;
    f = fopen("output.csv", "a+");
    char buffer[256]; 

    if (f == NULL) {
        printf("Error opening output file!\n");
        return 1;
    }
    

    // Create MPI_Datatype for Distance and Point structs
    create_mpi_point_type(&mpi_point_type);
    create_mpi_distance_type(&mpi_distance_type);
    
    double max_time_sum = 0.0;
   
    for (int iter = 0; iter < num_iterations; iter++) {
    
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
        MPI_Bcast(points, n, mpi_point_type, 0, MPI_COMM_WORLD);

        // Divide the work among processes
        int points_per_process = n / size;
        int start = rank * points_per_process;
        int end = (rank + 1) * points_per_process;
        if (rank == size - 1) end = n;  // Last process handles remaining points

        printf("Process %d: points_per_process = %d, start = %d, end = %d\n", rank, points_per_process, start, end);

        local_results = (Distance*)malloc(points_per_process * k * sizeof(Distance));
        printf("Process %d: Allocating memory for local_results[%d]\n", rank, points_per_process * k);

        start_time = MPI_Wtime();  // Start measuring parallel time

        // Each process finds k-nearest neighbors for its share of points
        find_k_nearest_neighbors(local_results, points, n, k, start, end);

        // Gather results
        if (rank == 0) {
            global_results = (Distance*)malloc(n * k * sizeof(Distance));
        }
        MPI_Gather(local_results, points_per_process * k, mpi_distance_type,
                global_results, points_per_process * k, mpi_distance_type,
                0, MPI_COMM_WORLD);
                
        //TODO
        /*if (rank == p - 1) {
            int remaining_points = n - (p - 1) * points_per_process;
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
            //printf("Parallel time taken: %f seconds\n", max_time);
            //fprintf(f, "Parallel time taken: %f seconds\n", max_time);
            //printf("\n%d Points - %d K neighbors - Average time over 1 iteration: %f seconds\n", n, k, max_time);
            //fprintf(f, "%d Points - %d K neighbors - Average time over 1 iteration: %f seconds\n", n, k, max_time);

            max_time_sum = max_time_sum + max_time;
            free(global_results);
        }
    
        free(local_results);
        free(points);

    }
    if(rank == 0){
        if (fgets(buffer, sizeof(buffer), f) == NULL) {
            fprintf(f, "Numero di processori (P),Numero di punti (N),Numero di vicini (K),Iterazioni,Tempo di esecuzione(s)\n");
        }
        printf("\n %d Num.Processor -%d Points - %d K neighbors - Average time over %d iterations: %f seconds\n", size, n, k, num_iterations, (max_time_sum/num_iterations));
        fprintf(f,"%d,%d,%d,%d,%f\n", size, n, k, num_iterations, (max_time_sum/num_iterations));
    }
    // Free the custom MPI_Datatypes
    MPI_Type_free(&mpi_point_type);
    MPI_Type_free(&mpi_distance_type);

    MPI_Finalize();
    fclose(f);
    return 0;
}

