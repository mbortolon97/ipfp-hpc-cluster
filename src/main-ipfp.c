#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "sparse_matrix.h"
#include "permutation.h"
#include "matrix_partition.h"
#include "submatrix.h"
#include "dense_distribution.h"
#include "log.h"

#define NUM_ITERATIONS 5


/*   100  40
    _________
140 | 50  20


*/




int main(int argc, char** argv) {
    if (argc != 5) {
        printf("Usage: distributedIPFP [aggregate_visit_matrix] [week_poi_marginals] [week_cbg_marginals] [output_dir]\n");
        exit(1);
    }
    srand(0);
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

    char log_filename[120];
    sprintf(log_filename, "log_%d.txt", world_rank);
    // load file
    FILE *in_file = fopen(log_filename, "w");

    // check file was opened correctly
    if (in_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    log_add_fp(in_file, 0);
    
    if(check_number_of_processes(&world_size, &world_rank)){
        MPI_Finalize();
        return 0; // the last process terminates if it won't be used
    }

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

    int hours;

    if (world_rank == 0) {
        // load matrices from files
        aggregate_visit_matrix = load_double_sparse_matrix(argv[1]);
        poi_marginals_matrix = load_double_sparse_matrix(argv[2]);
        cbg_marginals_matrix = load_double_dense_matrix(argv[3]);

        log_trace("aggregate_visit_matrix: %d %d\n", aggregate_visit_matrix.n_rows, aggregate_visit_matrix.n_cols);
        log_trace("poi_marginals_matrix: %d %d\n", poi_marginals_matrix.n_rows, poi_marginals_matrix.n_cols);
        log_trace("cbg_marginals_matrix: %d %d\n", cbg_marginals_matrix.n_rows, cbg_marginals_matrix.n_cols);

        // permute and get partitions
        permutation = create_sparse_matrix_random_permutation(aggregate_visit_matrix);
        permutated_aggregate_visit_matrix = permutate_double_sparse_matrix(permutation, aggregate_visit_matrix);
        partition = create_submatrix_partition(world_size, permutated_aggregate_visit_matrix.n_rows, permutated_aggregate_visit_matrix.n_cols);

        // send submatrices to other processes
        submatrix_to_elaborate = distribute_sparse_matrix(&partition, aggregate_visit_matrix); // TODO: replace with the permutated
        
        // clean_sparse_matrix(&aggregate_visit_matrix); /////////////////////////////////////////////////////////////////////////////////
        hours = poi_marginals_matrix.n_cols;
    } else {
        // receive submatrix from process_0
        submatrix_to_elaborate = wait_for_sparse_matrix();
    }

    log_trace("%d start_col: %d %d %d %d\n", world_rank, submatrix_to_elaborate.start_col, submatrix_to_elaborate.stop_col, submatrix_to_elaborate.start_row, submatrix_to_elaborate.stop_row);

    MPI_Bcast(&hours, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // for every hour compute IPFP
    for (int i = 0; i < hours; i++) {
        // vectors that will be used in IPFP (marginal sums of rows/cols)
        
        double_dense_matrix poi_marginals_at_hour_responsible;
        process_list col_process_list;
        double_dense_matrix cbg_marginals_at_hour_responsible;
        process_list row_process_list;

        if (world_rank == 0) {
            // get the column, apply the same permutation as before, broadcast
            double_dense_matrix poi_marginals_at_hour_not_perm = get_col_as_dense(poi_marginals_matrix, i);
            log_trace("%d poi_marginals_at_hour_not_perm:\n", world_rank);
            print_dense_matrix(poi_marginals_at_hour_not_perm);///////////////////////////////////////////////////////////////////
            double_dense_matrix cbg_marginals_at_hour_not_perm = get_row_from_dense(cbg_marginals_matrix, i);
            log_trace("%d cbg_marginals_at_hour_not_perm:\n", world_rank);
            print_dense_matrix(cbg_marginals_at_hour_not_perm);///////////////////////////////////////////////////////////////////
            double_dense_matrix cbg_marginals_at_hour = permutate_double_dense_matrix_along_rows(permutation, cbg_marginals_at_hour_not_perm);
            double_dense_matrix poi_marginals_at_hour = permutate_double_dense_matrix_along_columns(permutation, poi_marginals_at_hour_not_perm);
            // clean_double_dense_matrix(&poi_marginals_at_hour_not_perm); ////////////////////////////////////////////////////////////////////////////////////////////////
            // clean_double_dense_matrix(&cbg_marginals_at_hour_not_perm); ////////////////////////////////////////////////////////////////////////////////////////////////

            poi_marginals_at_hour_responsible = distribute_double_dense_matrix_using_column_partition(partition, poi_marginals_at_hour_not_perm, world_rank, MPI_COMM_WORLD); // TODO: replace with permutated
            clean_double_dense_matrix(&poi_marginals_at_hour);
            col_process_list = distribute_col_processes_list(partition, world_rank, MPI_COMM_WORLD);
            cbg_marginals_at_hour_responsible = distribute_double_dense_matrix_using_row_partition(partition, cbg_marginals_at_hour_not_perm, world_rank, MPI_COMM_WORLD); // TODO: replace with permutated
            clean_double_dense_matrix(&cbg_marginals_at_hour);
            row_process_list = distribute_row_processes_list(partition, world_rank, MPI_COMM_WORLD);
        } else {
            // receive broadcast
            if (row_responsible(submatrix_to_elaborate, world_rank)) {
                poi_marginals_at_hour_responsible = receive_double_dense_matrix(0, MPI_COMM_WORLD);
                row_process_list = receive_process_list(0, MPI_COMM_WORLD);
            }
            if (col_responsible(submatrix_to_elaborate, world_rank)) {
                cbg_marginals_at_hour_responsible = receive_double_dense_matrix(0, MPI_COMM_WORLD);
                col_process_list = receive_process_list(0, MPI_COMM_WORLD);
            }
        }
        // create a copy of the submatrix (one copy for each hour) 
        submatrix working_submatrix = clone_submatrix(submatrix_to_elaborate);

        log_trace("%d row_responsible %d\n", world_rank, working_submatrix.row_responsible);
        log_trace("%d col_responsible %d\n", world_rank, working_submatrix.col_responsible);

        int e = 0;

        for (int i = 1; i < NUM_ITERATIONS; i++) {
            log_trace("[%d] ITERATION %d\n",world_rank, i); ///////////////////////////////////////////////////////////////////
            util_print_submatrix(working_submatrix);////////////////////////////////////////////////////////////////////

            if (i % 2 == 1) {
                // last_w_axis_sum = sparse.sum(last_w, dim=0).to_dense()
                double_dense_matrix sum_result = sum_submatrix_along_cols(working_submatrix); // TODO: replace rows with cols without changing implementation
                log_trace("%d Sum result:\n", world_rank);
                print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                double_dense_matrix alfa_i;
                
                if (col_responsible(working_submatrix, world_rank)) {
                    log_trace("%d col_process_list.num_subprocesses %d\n", world_rank, col_process_list.num_subprocesses);
                    log_trace("%d Aggregate sum result (before):\n", world_rank);
                    print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                    aggregate_sum_results(sum_result, col_process_list, MPI_COMM_WORLD);
                    log_trace("%d Aggregate sum result:\n", world_rank);
                    print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                    // last_w_axis_sum[last_w_axis_sum < sys.float_info.epsilon] = 1.0
                    set_to_one_less_than_epsilon(sum_result);
                    log_trace("%d Set one less than epsilon:\n", world_rank);
                    print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                    // alfa_i = cbg_marginals_u / last_w_axis_sum
                    log_trace("%d Cbg marginal at hour responsible:\n", world_rank);
                    print_dense_matrix(cbg_marginals_at_hour_responsible);///////////////////////////////////////////////////////////////////
                    alfa_i = elementwise_division(cbg_marginals_at_hour_responsible, sum_result);
                    
                    log_trace("%d Elementwise division:\n", world_rank);
                    print_dense_matrix(alfa_i);///////////////////////////////////////////////////////////////////
                    
                    distribute_dense_matrix_to_processes(alfa_i, col_process_list, MPI_COMM_WORLD);
                } else {
                    send_sum_results(sum_result, working_submatrix.col_responsible, MPI_COMM_WORLD);
                    alfa_i = receive_double_dense_matrix(working_submatrix.col_responsible, MPI_COMM_WORLD);
                }
                if (world_rank==0) log_trace("[0] before multipl by alpha %lf\n",working_submatrix.elements[0].val); ///////////////////////////////////////////////////////////////////
                clean_double_dense_matrix(&sum_result);
                // new_w = sparse_dense_vector_mul(last_w, alfa_i)
                log_trace("alfa_i after distribution:\n");
                print_dense_matrix(alfa_i);///////////////////////////////////////////////////////////////////
                multiply_coefficient_by_cols(alfa_i, working_submatrix);
                printf("Multiply coefficient by cols:");
                util_print_submatrix(working_submatrix);////////////////////////////////////////////////////////////////////
                clean_double_dense_matrix(&alfa_i);
            } else {
                // last_w_axis_sum = torch.sparse.sum(last_w, dim=1).to_dense()
                double_dense_matrix sum_result = sum_submatrix_along_rows(working_submatrix); // TODO: replace cols with rows without changing implementation
                log_trace("%d Sum result:\n", world_rank);
                print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                double_dense_matrix alfa_i;
                if (row_responsible(working_submatrix, world_rank)) {
                    log_trace("%d row_process_list.num_subprocesses %d\n", world_rank, row_process_list.num_subprocesses);
                    log_trace("%d Aggregate sum result (before):\n", world_rank);
                    print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                    aggregate_sum_results(sum_result, row_process_list, MPI_COMM_WORLD);
                    log_trace("%d Aggregate sum result:\n", world_rank);
                    print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                    // last_w_axis_sum[last_w_axis_sum < sys.float_info.epsilon] = 1.0
                    set_to_one_less_than_epsilon(sum_result);
                    log_trace("%d Set one less than epsilon:\n", world_rank);
                    print_dense_matrix(sum_result);///////////////////////////////////////////////////////////////////
                    log_trace("%d POI marginal at hour responsible:\n", world_rank);
                    print_dense_matrix(poi_marginals_at_hour_responsible);///////////////////////////////////////////////////////////////////
                    // alfa_i = cbg_marginals_u / last_w_axis_sum
                    alfa_i = elementwise_division(poi_marginals_at_hour_responsible, sum_result);
                    log_trace("%d Elementwise division:\n", world_rank);
                    print_dense_matrix(alfa_i);///////////////////////////////////////////////////////////////////
                    distribute_dense_matrix_to_processes(alfa_i, row_process_list, MPI_COMM_WORLD);
                } else {
                    send_sum_results(sum_result, working_submatrix.row_responsible, MPI_COMM_WORLD);
                    alfa_i = receive_double_dense_matrix(working_submatrix.row_responsible, MPI_COMM_WORLD);
                }
                clean_double_dense_matrix(&sum_result);
                log_trace("alfa_i after distribution:\n");
                print_dense_matrix(alfa_i);///////////////////////////////////////////////////////////////////
                // new_w = sparse_dense_vector_mul(last_w, alfa_i)
                multiply_coefficient_by_rows(alfa_i, working_submatrix);
                printf("Multiply coefficient by rows:");
                util_print_submatrix(working_submatrix);////////////////////////////////////////////////////////////////////
                clean_double_dense_matrix(&alfa_i);
            }
        }

        if (row_responsible(submatrix_to_elaborate, world_rank)) {
            clean_double_dense_matrix(&poi_marginals_at_hour_responsible);
        }
        if (col_responsible(submatrix_to_elaborate, world_rank)) {
            clean_double_dense_matrix(&cbg_marginals_at_hour_responsible);
        }

        if (world_rank == 0) {
            // process_0 receives the submatrices, undo the permutations and store the result on a file
            double_sparse_matrix permutated_corrected_matrix = group_submatrices(working_submatrix, &aggregate_visit_matrix, partition.subp_cols*partition.subp_rows, partition.n_elements_biggest_partition);
            double_sparse_matrix corrected_matrix = inverse_double_sparse_matrix_permutation(permutation, permutated_corrected_matrix);
            clean_sparse_matrix(&permutated_corrected_matrix);
            char str[512];
            sprintf(str, "%s/save_result_%d.txt", argv[4], i);
            save_double_sparse_matrix(str, corrected_matrix);
            clean_sparse_matrix(&corrected_matrix);
        } else {
            // the other process send the submatrices to the main process
            send_submatrices(working_submatrix);
        }

        clean_submatrix(&working_submatrix);

        break;
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