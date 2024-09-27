#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include <stddef.h>
#include <string.h>

#define COORDINATOR 0

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

// Function to calculate Euclidean distance
float euclidean_distance(Point *point1, Point *point2) {
    float x = point1->x - point2->x;
    float y = point1->y - point2->y;
    float z = point1->z - point2->z;
    return sqrt(x * x + y * y + z * z);
}

// Function to create MPI datatype for Point structure
void create_mpi_point_type(MPI_Datatype *mpi_point_type) {
    int block_lengths[3] = {1, 1, 1};  // x, y, z are all floats
    MPI_Datatype types[3] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
    MPI_Aint offsets[3];

    // Calculate offsets of structure members
    offsets[0] = offsetof(Point, x);
    offsets[1] = offsetof(Point, y);
    offsets[2] = offsetof(Point, z);

    // Create MPI datatype for Point structure
    MPI_Type_create_struct(3, block_lengths, offsets, types, mpi_point_type);
    MPI_Type_commit(mpi_point_type);
}

// Function to create MPI datatype for Distance structure
void create_mpi_distance_type(MPI_Datatype *mpi_distance_type) {
    int block_lengths[2] = {1, 1};  // 1 float and 1 int
    MPI_Datatype types[2] = {MPI_FLOAT, MPI_INT};
    MPI_Aint offsets[2];

    // Calculate offsets of structure members
    offsets[0] = offsetof(Distance, neighbor_distance);
    offsets[1] = offsetof(Distance, neighbor_index);

    // Create MPI datatype for Distance structure
    MPI_Type_create_struct(2, block_lengths, offsets, types, mpi_distance_type);
    MPI_Type_commit(mpi_distance_type);
}

// Function to insert a neighbor into the K nearest neighbors list if it is closer
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

// Function to calculate distances and keep K nearest neighbors
void calculate_and_retain_k_nearest(Point *my_points, Point *received_points, 
                                    Distance *local_results, int points_per_process, 
                                    int K, int my_rank, int num_procs) {
    for (int i = 0; i < points_per_process; i++) {
        Distance *nearest_neighbors = &local_results[i * K];

        // Initialize the nearest neighbors array with max distance
        for (int k = 0; k < K; k++) {
            nearest_neighbors[k].neighbor_distance = INFINITY;
            nearest_neighbors[k].neighbor_index = -1;
        }

        // Calculate the distance of the current point to all received points
        for (int j = 0; j < points_per_process; j++) {
            int sender_rank = (my_rank - 1 + num_procs) % num_procs;
            int global_index = sender_rank * points_per_process + j;
            if (my_rank * points_per_process + i != global_index) {
                Distance new_neighbor;
                new_neighbor.neighbor_distance = euclidean_distance(&my_points[i], &received_points[j]);
                new_neighbor.neighbor_index = global_index;

                // Insert this new neighbor into the K nearest neighbors list if it is closer
                insert_distance_if_needed(nearest_neighbors, K, new_neighbor);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int my_rank, num_procs, N, K;
    Point *points = NULL;
    Point *my_points = NULL;
    Point *received_points = NULL;
    Distance *local_results = NULL;
    Distance *global_results = NULL;
    MPI_Datatype mpi_point_type, mpi_distance_type;
    int next, previous;
    double max_time;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (argc < 3 || argc > 4) {
        if (my_rank == COORDINATOR) {
            printf("Please provide a correct number of parameters\n");
            printf("Usage: <number_of_points> <number_of_neighbors> [-ev]\n");
            printf("Use [-ev] if you want to run the algorithm in evaluation mode\n");
        }
        MPI_Finalize();
        return 1;
    }

    N = atoi(argv[1]);
    K = atoi(argv[2]);

    next = (my_rank + 1) % (num_procs);
    previous = (my_rank - 1 + num_procs) % (num_procs);

    if (K >= N) {
        if (my_rank == COORDINATOR) {
            printf("K must be less than N.\n");
        }
        MPI_Finalize();
        return 1;
    }

    // Evaluation mode will run the code 5 times
    int evaluation_mode = (argc == 4 && strcmp(argv[3], "-ev") == 0) ? 1 : 0;
    int num_iterations = evaluation_mode ? 5 : 1;

    int points_per_process = N / num_procs;
    
    // Adjust for remainder if N is not exactly divisible by num_procs
    int remainder = N % num_procs;
    if (remainder != 0) {
        if (my_rank == COORDINATOR) {
            printf("Number of points must be divisible by the number of processes.\n");
        }
        MPI_Finalize();
        return 1;
    }

    FILE *f;
    f = fopen("output_ring.csv", "a+");
    char buffer[256]; 

    if (f == NULL) {
        printf("Error opening output file!\n");
        return 1;
    }

    // Create MPI datatypes for Point and Distance structures
    create_mpi_point_type(&mpi_point_type);
    create_mpi_distance_type(&mpi_distance_type);

    double max_time_sum = 0.0;
   
    for (int iter = 0; iter < num_iterations; iter++) {

        // Allocate memory for local points and K nearest neighbors
        my_points = (Point *)malloc(points_per_process * sizeof(Point));
        if (my_points == NULL) {
            printf("Process %d: Failed to allocate memory for my_points.\n", my_rank);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        local_results = (Distance *)malloc(points_per_process * K * sizeof(Distance));
        if (local_results == NULL) {
            printf("Process %d: Failed to allocate memory for local_results.\n", my_rank);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (my_rank == COORDINATOR) {
            points = (Point *)malloc(N * sizeof(Point));
            if (points == NULL) {
                printf("Coordinator: Failed to allocate memory for points.\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            srand(time(0));
            generate_points(points, N);
        }

        // Scatter the points to processes using custom MPI datatype
        MPI_Scatter(points, points_per_process, mpi_point_type,
                    my_points, points_per_process, mpi_point_type,
                    COORDINATOR, MPI_COMM_WORLD);

        // Initialize buffer for received points in the ring communication
        received_points = (Point *)malloc(points_per_process * sizeof(Point));
        if (received_points == NULL) {
            printf("Process %d: Failed to allocate memory for received_points.\n", my_rank);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        // Start time for measuring execution time
        double start_time = MPI_Wtime();

        for (int step = 0; step < num_procs; step++) {
            printf("Process %d: Starting computation step\n", my_rank);
            // Send my_points to the next process and receive points from the previous process using custom MPI datatype
            MPI_Sendrecv(my_points, points_per_process, mpi_point_type,
                        next, 0,
                        received_points, points_per_process, mpi_point_type,
                        previous, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Process %d: MPI_Sendrecv completed\n", my_rank);

            // Calculate distances and retain K nearest neighbors
            calculate_and_retain_k_nearest(my_points, received_points, local_results, points_per_process, K, my_rank, num_procs);
        }

        // Gather the results from all processes using custom MPI datatype
        if (my_rank == COORDINATOR) {
            global_results = (Distance *)malloc(N * K * sizeof(Distance));
            if (global_results == NULL) {
                printf("Coordinator: Failed to allocate memory for global_results.\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
        MPI_Gather(local_results, points_per_process * K, mpi_distance_type,
                global_results, points_per_process * K, mpi_distance_type,
                COORDINATOR, MPI_COMM_WORLD);

        // End time for measuring execution time
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        max_time_sum = max_time_sum + max_time;

        // Only the coordinator process will print the results
        if (my_rank == COORDINATOR) {
            /*for (int i = 0; i < N; i++) {
                printf("Point %d's %d nearest neighbors:\n", i, K);
                for (int m = 0; m < K; m++) {
                    printf("Neighbor %d: Point %d (Distance: %f)\n", m + 1, global_results[i * K + m].neighbor_index, global_results[i * K + m].neighbor_distance);
                }
            }*/
            printf("\n %d Num.Processor - %d Points - %d K neighbors - Time: %f seconds\n", num_procs, N, K, max_time);
        }

        // Free memory and MPI datatypes
        if (my_rank == COORDINATOR) {
            free(points);
            free(global_results);
        }
        free(my_points);
        free(received_points);
        free(local_results);
    }
    if(my_rank == 0){
        if (fgets(buffer, sizeof(buffer), f) == NULL) {
            fprintf(f, "Numero di processori (P),Numero di punti (N),Numero di vicini (K),Iterazioni,Tempo di esecuzione(s)\n");
        }
        printf("\n %d Num.Processor -%d Points - %d K neighbors - Average time over %d iterations: %f seconds\n", num_procs, N, K, num_iterations, (max_time_sum/num_iterations));
        fprintf(f,"%d,%d,%d,%d,%f\n", num_procs, N, K, num_iterations, (max_time_sum/num_iterations));
    }
    MPI_Type_free(&mpi_point_type);
    MPI_Type_free(&mpi_distance_type);

    MPI_Finalize();
    return 0;
}
