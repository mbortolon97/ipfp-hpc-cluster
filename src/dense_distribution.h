#ifndef __DENSE_DISTRIBUTION_H
#define __DENSE_DISTRIBUTION_H 1

#include "dense_matrix.h"
#include <mpi.h>
#include "matrix_partition.h"

double_dense_matrix distribute_double_dense_matrix_using_column_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm);
double_dense_matrix distribute_double_dense_matrix_using_row_partition(const submatrix_partition partition, const double_dense_matrix matrix, int world_rank, MPI_Comm comm);
double_dense_matrix receive_dense_matrix(int source, MPI_Comm comm);

#endif