#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "input.h"
#include "submatrix.h"
#include "matrix_partition.h"



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


submatrix distribute_sparse_matrix(submatrix_partition partition, double_sparse_matrix matirx){
    int n_processes = partition.subp_rows * partition.subp_cols;
    int elements_assigned_to_process_counter[n_processes] = {0};
    MPI_Datatype mpi_tuple;
    submatrix process_0_submatrix;

    for (int i=0; i<matrix.n_elements; i++){
        int found = 0;
        for (int p=0; p<n_processes && !found; p++){
            // if process P must handle element I
            if (check_if_inside_submatrix(partition.assignments[p], matrix.rows[i], matrix.cols[i])){
                found=1;
                // increase the counter of the elements handled by P
                elements_assigned_to_process_counter[p]++;
            }
        }
    }

    // MPI datatype
    build_mpi_tuple(&mpi_tuple);

    // send infos to other processes
    for (int p=1; p<n_processes && !found; p++){
        
        // send some informations to the specific process (so that he can allocate the correct space in memory)
        int infos[7];
        infos[0] = partition.assignments[p].start_row;
        infos[1] = partition.assignments[p].stop_row;
        infos[2] = partition.assignments[p].start_col;
        infos[3] = partition.assignments[p].stop_col;
        infos[4] = elements_assigned_to_process_counter[p];
        infos[5] = partition.assignments[p].col_responsible;
        infos[6] = partition.assignments[p].row_responsible;
  
        // create a countiguous block of memory to send the data
        struct mpi_sparse_matrix_element* mpi_data = malloc(sizeof(mpi_sparse_matrix_element)*elements_assigned_to_process_counter[p]);
        int count = 0;
        for (int i=0; i<matrix.n_elements; i++){
            if (check_if_inside_submatrix(partition.assignments[p], matrix.rows[i], matrix.cols[i])){
                mpi_data[count].row = matrix.rows[i];
                mpi_data[count].col = matrix.cols[i];
                mpi_data[count].val = matrix.values[i];
                count++;
            }
        }
        if (p==0)
            process_0_submatrix = create_submatrix(); //// TO BE DONE
        else{
            MPI_Send( infos , 7 , MPI_INT , p , 0 , MPI_COMM_WORLD);
            MPI_Send( mpi_data , elements_assigned_to_process_counter[p] , mpi_tuple , p , 0 , MPI_COMM_WORLD);
        }
        free(mpi_sparse_matrix_element);
    }    

    return process_0_submatrix;
}

wait_for_sparse_matrix