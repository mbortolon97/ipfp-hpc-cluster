#ifndef __SUBMATRIX_H
#define __SUBMATRIX_H 1

#include "sparse_matrix.h"
#include "matrix_partition.h"

struct mpi_matrix_element {
    int row, col;
    double val;
};

typedef struct submatrix_queue_struct{
    struct submatrix_queue_struct* next;
    struct mpi_matrix_element* element;
} submatrix_queue;

typedef struct submatrix_struct
{
        int start_row;
        int stop_row;
        int start_col;
        int stop_col;
        
        int n_elements;
        struct mpi_matrix_element* elements;

        int col_responsible;
        int row_responsible;

        submatrix_queue** row_queue;
        submatrix_queue** col_queue;
} submatrix;


// get submatrix through MPI
submatrix distribute_sparse_matrix(submatrix_partition partition, double_sparse_matrix matirx);
submatrix wait_for_sparse_matrix();

// operations on submatrices
double_dense_matrix sum_submatrix_along_rows(const submatrix submatrix);
double_dense_matrix sum_submatrix_along_cols(const submatrix submatrix);

void multiply_coefficient_by_cols(const double_dense_matrix alfa_i, submatrix working_submatrix);
void multiply_coefficient_by_rows(const double_dense_matrix alfa_i, submatrix working_submatrix);

#endif