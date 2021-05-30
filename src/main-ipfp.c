#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sparse_matrix.h"
#include "permutation.h"
#include "matrix_partition.h"
#include "submatrix.h"
#include "dense_distribution.h"

#define NUM_ITERATIONS 100

/*
 * TODOS:
 * - if n_processes is a prime number >= 5 ---> test (we should kill the last process which is unused and decrease world size)
 * - if a submatrix is empty, does it still work?
 * - line 72: permutated_aggregate_visit_matrix ?
 * - remove submatrix_queue attribute from submatrix
 */


int main(int argc, char** argv) {
    if (argc != 4) {
        printf("Usage: distributedIPFP [aggregate_visit_matrix] [week_poi_marginals] [week_cbg_marginals]\n");
        exit(1);
    }
    srand(time(NULL));
    // Input ( read from the file - file format? )
    // Shuffle rows and columns (create map between origin and destination)
    // Compute the sub matrix divisions (greedy algorithm)

    // Assign each sub-matrix to a node
    // Assign which node is responsible for each row and column
    // Communicate this information to the various nodes
    // Send each submatrix to the appropriate node using the map of the shuffle
    
    // IPFP cycle (initially for 100 iterations):
    // - Compute the sum for the row/column inside the submatrix
    // - Send the result to the row master

	MPI_Init(&argc, &argv);
	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // input data read from files.txt
    double_sparse_matrix aggregate_visit_matrix;
    double_sparse_matrix poi_marginals_matrix;
    double_dense_matrix cbg_marginals_matrix;

    // aggregate_visit_matrix  after beign permutated
    double_sparse_matrix permutated_aggregate_visit_matrix;

    // the permutation used
    sparse_matrix_permutation permutation;

    // restructuring sparse_matrix_permutation to send it through MPI
    submatrix_partition partition; 
    
    // submatrix where the processes will work on
    submatrix submatrix_to_elaborate;

    if (world_rank == 0) {
        // load matrices from files
        aggregate_visit_matrix = load_double_sparse_matrix(argv[1]);
        poi_marginals_matrix = load_double_sparse_matrix(argv[2]);
        cbg_marginals_matrix = load_double_dense_matrix(argv[3]);

        // permute and get partitions
        permutation = create_sparse_matrix_random_permutation(aggregate_visit_matrix);
        permutated_aggregate_visit_matrix = permutate_double_sparse_matrix(permutation, aggregate_visit_matrix);
        partition = create_submatrix_partition(world_size, aggregate_visit_matrix.n_rows, aggregate_visit_matrix.n_cols);

        printf("Sending data...\n");
        // send submatrices to other processes
        submatrix_to_elaborate = distribute_sparse_matrix(partition, aggregate_visit_matrix);
        printf("Data send\n");
    } else {
        printf("Receiving data... %d\n", world_rank);
        // receive submatrix from process_0
        submatrix_to_elaborate = wait_for_sparse_matrix();
        printf("Data received %d\n", world_rank);
    }

    printf("Initialization success %d\n", world_rank);

    // for every hour compute IPFP
    for (int i = 0; i < poi_marginals_matrix.n_cols; i++) {
        // vectors that will be used in IPFP (marginal sums of rows/cols)

        double_dense_matrix poi_marginals_at_hour_responsible;
        process_list col_process_list;
        double_dense_matrix cbg_marginals_at_hour_responsible;
        process_list row_process_list;

        if (world_rank == 0) {
            // get the column, apply the same permutation as before, broadcast
            double_dense_matrix poi_marginals_at_hour_not_perm = get_col_as_dense(poi_marginals_matrix, i);
            printf("Col as dense\n");
            double_dense_matrix cbg_marginals_at_hour_not_perm = get_row_from_dense(cbg_marginals_matrix, i);
            printf("Row as dense\n");
            double_dense_matrix cbg_marginals_at_hour = permutate_double_dense_matrix_along_rows(permutation, cbg_marginals_at_hour_not_perm);
            printf("Permutate double dense matrix along rows");
            double_dense_matrix poi_marginals_at_hour = permutate_double_dense_matrix_along_columns(permutation, poi_marginals_at_hour_not_perm);
            printf("Permutate double dense matrix along columns");
            clean_double_dense_matrix(&poi_marginals_at_hour_not_perm);
            clean_double_dense_matrix(&cbg_marginals_at_hour_not_perm);

            poi_marginals_at_hour_responsible = distribute_double_dense_matrix_using_column_partition(partition, poi_marginals_at_hour, world_rank, MPI_COMM_WORLD);
            clean_double_dense_matrix(&poi_marginals_at_hour);
            col_process_list = distribute_col_processes_list(partition, world_rank, MPI_COMM_WORLD);
            cbg_marginals_at_hour_responsible = distribute_double_dense_matrix_using_row_partition(partition, cbg_marginals_at_hour, world_rank, MPI_COMM_WORLD);
            clean_double_dense_matrix(&cbg_marginals_at_hour);
            row_process_list = distribute_row_processes_list(partition, world_rank, MPI_COMM_WORLD);
            
        } else {
            // receive broadcast
            if (col_responsible(submatrix_to_elaborate, world_rank) == 0) {
                poi_marginals_at_hour_responsible = receive_dense_matrix(0, MPI_COMM_WORLD);
                col_process_list = receive_process_list(0, MPI_COMM_WORLD);
            }
            if (row_responsible(submatrix_to_elaborate, world_rank) == 0) {
                cbg_marginals_at_hour_responsible = receive_dense_matrix(0, MPI_COMM_WORLD);
                row_process_list = receive_process_list(0, MPI_COMM_WORLD);
            }
        }
        // create a copy of the submatrix (one copy for each hour) 
        submatrix working_submatrix = clone_submatrix(submatrix_to_elaborate);

        printf("Complete success");

        for (int i = 0; i < NUM_ITERATIONS; i++) {
            if (i % 2 == 1) {
                // last_w_axis_sum = sparse.sum(last_w, dim=0).to_dense()
                double_dense_matrix sum_result = sum_submatrix_along_rows(working_submatrix);
                double_dense_matrix alfa_i;
                if (row_responsible(working_submatrix, world_rank)) {
                    aggregate_sum_results(sum_result, row_process_list, MPI_COMM_WORLD);
                    // last_w_axis_sum[last_w_axis_sum < sys.float_info.epsilon] = 1.0
                    set_to_one_less_than_epsilon(sum_result);
                    // alfa_i = cbg_marginals_u / last_w_axis_sum
                    alfa_i = elementwise_division(cbg_marginals_at_hour_responsible, sum_result);
                    distribute_dense_matrix_to_processes(alfa_i, row_process_list, MPI_COMM_WORLD);
                } else {
                    send_sum_results(sum_result, working_submatrix.row_responsible, MPI_COMM_WORLD);
                    alfa_i = receive_dense_matrix(working_submatrix.row_responsible, MPI_COMM_WORLD);
                }
                clean_double_dense_matrix(&sum_result);
                // new_w = sparse_dense_vector_mul(last_w, alfa_i)
                multiply_coefficient_by_cols(alfa_i, working_submatrix);
                clean_double_dense_matrix(&alfa_i);
            } else {
                // last_w_axis_sum = torch.sparse.sum(last_w, dim=1).to_dense()
                double_dense_matrix sum_result = sum_submatrix_along_cols(working_submatrix);
                double_dense_matrix alfa_i;
                if (col_responsible(working_submatrix, world_rank)) {
                    aggregate_sum_results(sum_result, col_process_list, MPI_COMM_WORLD);
                    // last_w_axis_sum[last_w_axis_sum < sys.float_info.epsilon] = 1.0
                    set_to_one_less_than_epsilon(sum_result);
                    // alfa_i = cbg_marginals_u / last_w_axis_sum
                    alfa_i = elementwise_division(poi_marginals_at_hour_responsible, sum_result);
                    distribute_dense_matrix_to_processes(alfa_i, col_process_list, MPI_COMM_WORLD);
                } else {
                    send_sum_results(sum_result, working_submatrix.col_responsible, MPI_COMM_WORLD);
                    alfa_i = receive_dense_matrix(working_submatrix.col_responsible, MPI_COMM_WORLD);
                }
                clean_double_dense_matrix(&sum_result);
                // new_w = sparse_dense_vector_mul(last_w, alfa_i)
                multiply_coefficient_by_rows(alfa_i, working_submatrix);
                clean_double_dense_matrix(&alfa_i);
            }
        }
        if (col_responsible(submatrix_to_elaborate, world_rank)) {
            clean_double_dense_matrix(&poi_marginals_at_hour_responsible);
        }
        if (row_responsible(submatrix_to_elaborate, world_rank)) {
            clean_double_dense_matrix(&cbg_marginals_at_hour_responsible);
        }
        /*
        if (world_rank == 0) {
            // process_0 receives the submatrices, undo the permutations and store the result on a file
            double_sparse_matrix permutated_corrected_matrix = group_submatrices(working_submatrix, &aggregate_visit_matrix, partition.subp_cols*partition.subp_rows, partition.n_elements_biggest_partition);
            double_sparse_matrix corrected_matrix = inverse_double_sparse_matrix_permutation(permutation, permutated_corrected_matrix);
            clean_sparse_matrix(&permutated_corrected_matrix);
            char str[512];
            sprintf(str, "save_result_%d.txt", i);
            save_double_sparse_matrix(str, corrected_matrix);
            clean_sparse_matrix(&corrected_matrix);
            
        } else {
            // the other process send the submatrices to the main process
            send_submatrices(working_submatrix); 
        }
        */

        clean_submatrix(&working_submatrix);
    }
    
    if(world_rank == 0) {
        clean_sparse_matrix_permutation(&permutation);
        clean_sparse_matrix(&aggregate_visit_matrix);
        clean_sparse_matrix(&poi_marginals_matrix);
        clean_double_dense_matrix(&cbg_marginals_matrix);
    }
	
    MPI_Finalize();
    return 1;
}