#include "dense_distribution.h"
#include <assert.h>
#include <string.h>
#include "message_tag.h"

#include <stdlib.h>

double_dense_matrix distribute_double_dense_matrix_using_column_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm) {
    double_dense_matrix result_matrix;
    assert(matrix.n_rows == 1);

    int dimension[2];
    dimension[0] = matrix.n_rows;
    
    int i;
    for (int i = 0; i < partition.subp_cols; i++) {
        dimension[1] = partition.col_master[i].stop_col - partition.col_master[i].start_col;
        if (partition.col_master[i].col_master_process_id != world_rank) {
            MPI_Send(&dimension, 2, MPI_INT, partition.col_master[i].col_master_process_id, SEND_DENSE_MATRIX_DIMENSION, comm);
            MPI_Send(&(matrix.matrix[partition.col_master[i].start_col]), dimension[1], MPI_DOUBLE, partition.col_master[i].col_master_process_id, SEND_DENSE_MATRIX_CONTENT, comm);
        } else {
            result_matrix = create_double_dense_matrix(dimension[0], dimension[1]);
            memcpy(result_matrix.matrix, &(matrix.matrix[partition.col_master[i].start_col]), dimension[1] * sizeof(double));
        }
    }
    return result_matrix;
}

double_dense_matrix distribute_double_dense_matrix_using_row_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm) {
    double_dense_matrix result_matrix;
    assert(matrix.n_cols == 1);

    int dimension[2];
    dimension[1] = matrix.n_cols;

    int i;
    for (int i = 0; i < partition.subp_rows; i++) {
        dimension[0] = partition.row_master[i].stop_row - partition.row_master[i].start_row;
        if (partition.row_master[i].row_master_process_id != world_rank) {
            MPI_Send(&dimension, 2, MPI_INT, partition.row_master[i].row_master_process_id, SEND_DENSE_MATRIX_DIMENSION, comm);
            MPI_Send(&(matrix.matrix[partition.row_master[i].start_row]), dimension[0], MPI_DOUBLE, partition.row_master[i].row_master_process_id, SEND_DENSE_MATRIX_CONTENT, comm);
        } else {
            result_matrix = create_double_dense_matrix(dimension[0], dimension[1]);
            memcpy(result_matrix.matrix, &(matrix.matrix[partition.row_master[i].stop_row]), dimension[0] * sizeof(double));
        }
    }
}

process_list distribute_row_processes_list(const submatrix_partition partition, int world_rank, MPI_Comm comm) {
    process_list list;
    int i;
    for (i = 0; i < partition.subp_rows; i++) {
        if (partition.row_master[i].row_master_process_id == world_rank) {
            list.num_subprocesses = partition.row_master[i].num_subprocesses;
            list.processes_id = malloc(sizeof(int) * list.num_subprocesses);
            memcpy(list.processes_id, partition.row_master[i].row_processes_id, sizeof(int) * list.num_subprocesses);
        } else {
            MPI_Send(&(partition.row_master[i].num_subprocesses), 1, MPI_INT, partition.row_master[i].row_master_process_id, SEND_PROCESS_LIST_N_PROC, comm);
            MPI_Send(&(partition.row_master[i].row_processes_id), partition.row_master[i].num_subprocesses, MPI_INT, partition.row_master[i].row_master_process_id, SEND_PROCESS_LIST_PROCESSES, comm);
        }
    }
    return list;
}

process_list distribute_col_processes_list(const submatrix_partition partition, int world_rank, MPI_Comm comm) {
    process_list list;
    int i;
    for (i = 0; i < partition.subp_cols; i++) {
        if (partition.col_master[i].col_master_process_id == world_rank) {
            list.num_subprocesses = partition.col_master[i].num_subprocesses;
            list.processes_id = malloc(sizeof(int) * list.num_subprocesses);
            memcpy(list.processes_id, partition.col_master[i].col_processes_id, sizeof(int) * list.num_subprocesses);
        } else {
            MPI_Send(&(partition.col_master[i].num_subprocesses), 1, MPI_INT, partition.col_master[i].col_master_process_id, SEND_PROCESS_LIST_N_PROC, comm);
            MPI_Send(&(partition.col_master[i].col_processes_id), partition.col_master[i].num_subprocesses, MPI_INT, partition.col_master[i].col_master_process_id, SEND_PROCESS_LIST_PROCESSES, comm);
        }
    }
    return list;
}

process_list receive_process_list(int source, MPI_Comm comm) {
    MPI_Status status;
    process_list list;
    int n_subprocess;
    MPI_Recv(&(list.num_subprocesses), 1, MPI_INT, source, SEND_PROCESS_LIST_N_PROC, comm, &status);
    list.processes_id = malloc(sizeof(int) * list.num_subprocesses);
    MPI_Recv(list.processes_id, partition.col_master[i].num_subprocesses, MPI_INT, partition.col_master[i].col_master_process_id, SEND_PROCESS_LIST_PROCESSES, comm);
    return list;
}

void clean_process_list(process_list* list) {
    free(list->processes_id);
    list->processes_id = NULL;
}

double_dense_matrix receive_dense_matrix(int source, MPI_Comm comm) {
    MPI_Status status;
    int dimension[2];
    MPI_Recv(dimension, 2, MPI_INT, source, SEND_DENSE_MATRIX_DIMENSION, comm, &status);
    double_dense_matrix result_matrix = create_double_dense_matrix(dimension[0], dimension[1]);
    MPI_Recv(result_matrix.matrix, 2, MPI_INT, source, SEND_DENSE_MATRIX_CONTENT, comm, &status);
    return result_matrix;
}

void aggregate_sum_results(double_dense_matrix sum_result, const process_list list, MPI_Comm comm) {
    MPI_Status status;
    double_dense_matrix receiving_buffer = create_double_dense_matrix(sum_result.n_rows, sum_result.n_cols);
    int i, j;
    for (i = 0; i < list.num_subprocesses; i++) {
        MPI_Recv(receiving_buffer.matrix, sum_result.n_rows * sum_result.n_cols, MPI_INT, MPI_ANY_SOURCE, SEND_SUM_RESULTS, comm, &status);

        for (j = 0; j < result_matrix.matrix; j++) {
            sum_result.matrix[j] += receiving_buffer.matrix[j];
        }
    }
    clean_double_dense_matrix(&receiving_buffer);
}

void send_sum_results(double_dense_matrix sum_result, int dest, MPI_Comm comm) {
    MPI_Send(sum_result.matrix, sum_result.n_rows * sum_result.n_cols, MPI_INT, dest, SEND_SUM_RESULTS, comm);
}

void distribute_dense_matrix_to_processes(double_dense_matrix matrix, process_list list, MPI_Comm comm) {
    int dimension[2];
    dimension[0] = matrix.n_rows;
    dimension[1] = matrix.n_cols;
    int i;
    for (i = 0; i < list.num_subprocesses; i++) {
        MPI_Send(&dimension, 2, MPI_INT, list.processes_id[i], SEND_DENSE_MATRIX_DIMENSION, comm);
        MPI_Send(matrix.matrix, dimension[0] * dimension[1], MPI_DOUBLE, list.processes_id[i], SEND_DENSE_MATRIX_CONTENT, comm);
    }
}