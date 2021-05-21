#include <stdio.h>
#include "input.h"
#include "mpiutils.h"

int main()
{
    MPI_Init(NULL, NULL);
    setbuf(stdout, NULL);

    struct sparse_matrix *m;
    int elements_per_process[300];
    struct sparse_matrix_element *mpi_matrix;
    MPI_Datatype mpi_tuple;

    m = load_sparse_matrix();
    printf("rows %d,  elements %d\n", m->n_rows, m->n_elements);

    build_mpi_tuple(&mpi_tuple);

    mpi_matrix = split_in_submatrices(m, 100, elements_per_process);
    printf("ended\n");

    // scatter

    //queue matrix to make sum fast

    // how to send messages

    // re-build original matrix

    MPI_Finalize();

    return 0;
}
