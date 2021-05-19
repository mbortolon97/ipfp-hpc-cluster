#include <mpi.h>
#include "input.h"

#ifndef MPIUTILS
#define MPIUTILS

void split_in_submatrices(struct sparse_matrix *, int n_processes);

#endif