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
} submatrix;

// returns 0 if you are a col/row master
int col_responsible(const submatrix submatrix, const int world_rank);
int row_responsible(const submatrix submatrix, const int world_rank);

// get submatrix through MPI
submatrix distribute_sparse_matrix(submatrix_partition partition, double_sparse_matrix matirx);
submatrix wait_for_sparse_matrix();

submatrix clone_submatrix(const submatrix original_submatrix);

// send submatrices through MPI
double_sparse_matrix group_submatrices(submatrix working_submatrix, double_sparse_matrix* original_matrix, const int n_processes, const int n_elements_biggest_sumbatrix);
void send_submatrices(submatrix working_submatrix); 



// operations on submatrices
double_dense_matrix sum_submatrix_along_rows(const submatrix submatrix);
double_dense_matrix sum_submatrix_along_cols(const submatrix submatrix);

void multiply_coefficient_by_cols(const double_dense_matrix alfa_i, submatrix working_submatrix);
void multiply_coefficient_by_rows(const double_dense_matrix alfa_i, submatrix working_submatrix);

#endif