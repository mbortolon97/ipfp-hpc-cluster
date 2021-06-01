#ifndef __SUBMATRIX_H
#define __SUBMATRIX_H 1

#include "sparse_matrix.h"
#include "matrix_partition.h"
#include <stdlib.h>

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

/**
 * This function checks if the process is the col master
 **/ 
static inline bool col_responsible(const submatrix submatrix, const int world_rank) {
    return submatrix.col_responsible == world_rank;
}

/**
 * This function checks if the process is the row master
 **/ 
static inline bool row_responsible(const submatrix submatrix, const int world_rank) {
    return submatrix.row_responsible == world_rank;
}


// distribute submatrix through MPI

/**
 * This function distributes the aggregate matrix across the subprocesses, returns its own submatrix
 **/ 
submatrix distribute_sparse_matrix(submatrix_partition *partition, double_sparse_matrix matirx);

/**
 * This function returns the submatrix received from process 0
 **/ 
submatrix wait_for_sparse_matrix();


/**
 * This function makes a copy of the submatrix (one copy for each hour)
 **/ 
submatrix clone_submatrix(const submatrix original_submatrix);


// regroup submatrices through MPI

/**
 * This function gather the submatrices across the processes to have a final matrix in process 0
 **/ 
double_sparse_matrix group_submatrices(submatrix working_submatrix, double_sparse_matrix* original_matrix, const int n_processes, const int n_elements_biggest_sumbatrix);

/**
 * This function send the submatrix to process 0
 **/ 
void send_submatrices(submatrix working_submatrix); 


// operations on submatrices

/**
 * This function computes the partial sum along the rows of a matrix, returns a column vector
 **/ 
double_dense_matrix sum_submatrix_along_rows(const submatrix submatrix);

/**
 * This function computes the partial sum along the cols of a matrix, returns a row vector
 **/ 
double_dense_matrix sum_submatrix_along_cols(const submatrix submatrix);

/**
 * This function does the inplace-broadcasted-elementwise multiplication between the alphas shape(1,ncol) and the submatrix shape(nrow,ncol)
 **/ 
void multiply_coefficient_by_cols(const double_dense_matrix alfa_i, submatrix working_submatrix);

/**
 * This function does the inplace-broadcasted-elementwise multiplication between the alphas shape(nrow,1) and the submatrix shape(nrow,ncol)
 **/ 
void multiply_coefficient_by_rows(const double_dense_matrix alfa_i, submatrix working_submatrix);

/**
 * This function frees the memory allocated for a submatrix
 **/ 
void clean_submatrix(submatrix *submatrix);


// utils

/**
 * Just a function that writes in stdout the content of a submatrix
 **/ 
void util_print_submatrix(const submatrix submatrix);

#endif