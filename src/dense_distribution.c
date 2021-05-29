#include "dense_distribution.h"
#include <assert.h>
#include <string.h>
#include "message_tag.h"

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

double_dense_matrix receive_dense_matrix(int source, MPI_Comm comm) {
    MPI_Status status;
    int dimension[2];
    MPI_Recv(dimension, 2, MPI_INT, source, SEND_DENSE_MATRIX_DIMENSION, comm, &status);
    double_dense_matrix result_matrix = create_double_dense_matrix(dimension[0], dimension[1]);
    MPI_Recv(result_matrix.matrix, 2, MPI_INT, source, SEND_DENSE_MATRIX_CONTENT, comm, &status);
    return result_matrix;
}