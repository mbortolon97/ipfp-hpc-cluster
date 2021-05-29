#include "matrix_partition.h"
#include <mpi.h>

#ifndef __SUBMATRIX_H
#define __SUBMATRIX_H 1


struct mpi_sparse_matrix_element {
    int row, col;
    double val;
};

typedef struct submatrix_queue_struct{
    submatrix_queue* next;
    struct mpi_sparse_matrix_element* element;
} submatrix_queue;

typedef struct submatrix_struct
{
        int start_row;
        int stop_row;
        int start_col;
        int stop_col;
        
        int n_elements;
        struct mpi_sparse_matrix_element* elements;

        int col_responsible;
        int row_responsible;

        submatrix_queue** row_queue;
        submatrix_queue** col_queue;
} submatrix;



void build_mpi_tuple(MPI_Datatype *mpi_tuple)
{
    struct mpi_sparse_matrix_element element;
    // array of structure member sizes
    int blocklengths[3] = {1, 1, 1};

    // structure member types
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_FLOAT};

    // offset of structure members
    MPI_Aint offsets[3];
    MPI_Aint points[3];
    MPI_Get_address(&element.row, points);
    MPI_Get_address(&element.col, points + 1);
    MPI_Get_address(&element.val, points + 2);
    offsets[0] = 0;
    offsets[1] = points[1] - points[0];
    offsets[2] = points[2] - points[0];

    // create mpi struct
    MPI_Type_create_struct(3, blocklengths, offsets, types, mpi_tuple);
    MPI_Type_commit(mpi_tuple);
}

// get submatrix through MPI
submatrix distribute_sparse_matrix(submatrix_partition partition, double_sparse_matrix matirx);
submatrix wait_for_sparse_matrix();

// operations on submatrices




#endif