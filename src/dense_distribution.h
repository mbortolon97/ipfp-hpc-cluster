#ifndef __DENSE_DISTRIBUTION_H
#define __DENSE_DISTRIBUTION_H 1

#include "dense_matrix.h"
#include <mpi.h>
#include "matrix_partition.h"

typedef struct process_list_struct {
    int num_subprocesses;
    int* processes_id;
} process_list;

// functions that handle communication of marginals across processes

/**
 * This function returns a new dense matrix containing the marginals of its own column partition and sends the other marginals to the other master processes
 **/ 
double_dense_matrix distribute_double_dense_matrix_using_column_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm);

/**
 * This function returns a new dense matrix containing the marginals of its own row partition and sends the other marginals to the other master processes
 **/ 
double_dense_matrix distribute_double_dense_matrix_using_row_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm);

/**
 * This function returns a new dense matrix containing the marginals, which is received from process 0
 **/ 
double_dense_matrix receive_double_dense_matrix(int source, MPI_Comm comm);


// functions that handle communication of list of subprocesses across processes

/**
 * This function returns a list of the subprocesses handled by process 0 and sends the list to the other row masters
 **/ 
process_list distribute_row_processes_list(const submatrix_partition partition, int world_rank, MPI_Comm comm);

/**
 * This function returns a list of the subprocesses handled by process 0 and sends the list to the other col masters
 **/ 
process_list distribute_col_processes_list(const submatrix_partition partition, int world_rank, MPI_Comm comm);

/**
 * This function returns a list of the subprocesses, received from process 0
 **/ 
process_list receive_process_list(int source, MPI_Comm comm);

/**
 * free the memory allocated for a list
 **/ 
void clean_process_list(process_list* list);


// functions that handle communication of partial sums across processes

/**
 * This function sends the partial sum to the row/col master
 **/ 
void send_sum_results(double_dense_matrix sum_result, int dest, MPI_Comm comm);

/**
 * This function returns the final sum of the partial sums received by the subprocesses
 **/ 
void aggregate_sum_results(double_dense_matrix sum_result, const process_list list, MPI_Comm comm);


// functions that handle communication of alpha along rows/columns

/**
 * This function sends the alphas to each subprocess in the list of a master
 **/
void distribute_dense_matrix_to_processes(double_dense_matrix matrix, process_list list, MPI_Comm comm);

#endif