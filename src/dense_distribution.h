#ifndef __DENSE_DISTRIBUTION_H
#define __DENSE_DISTRIBUTION_H 1

#include "dense_matrix.h"
#include <mpi.h>
#include "matrix_partition.h"

typedef struct process_list_struct {
    int num_subprocesses;
    int* processes_id;
} process_list;

double_dense_matrix distribute_double_dense_matrix_using_column_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm);
double_dense_matrix distribute_double_dense_matrix_using_row_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm);
double_dense_matrix receive_double_dense_matrix(int source, MPI_Comm comm);


process_list distribute_row_processes_list(const submatrix_partition partition, int world_rank, MPI_Comm comm);
process_list distribute_col_processes_list(const submatrix_partition partition, int world_rank, MPI_Comm comm);
process_list receive_process_list(int source, MPI_Comm comm);
void clean_process_list(process_list* list);

void send_sum_results(double_dense_matrix sum_result, int dest, MPI_Comm comm);
void aggregate_sum_results(double_dense_matrix sum_result, const process_list list, MPI_Comm comm);

void distribute_dense_matrix_to_processes(double_dense_matrix matrix, process_list list, MPI_Comm comm);

#endif