#include <mpi.h>
#include "input.h"

#ifndef MPIUTILS
#define MPIUTILS

struct sparse_matrix_element *split_in_submatrices(struct sparse_matrix *matrix, int n_processes, int *elements_per_process);

void build_mpi_tuple(MPI_Datatype *mpi_tuple);

#endif