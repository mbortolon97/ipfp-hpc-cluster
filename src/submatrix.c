#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "submatrix.h"
#include "matrix_partition.h"

submatrix create_empty_submatrix(int* infos);

void build_mpi_tuple(MPI_Datatype *mpi_tuple)
{
    struct mpi_matrix_element element;
    // array of structure member sizes
    int blocklengths[3] = {1, 1, 1};

    // structure member types
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_DOUBLE};

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


submatrix distribute_sparse_matrix(submatrix_partition *partition, double_sparse_matrix matrix){
    int n_processes = partition->subp_rows * partition->subp_cols;
    int elements_assigned_to_process_counter[n_processes];
    MPI_Datatype mpi_tuple;
    submatrix process_0_submatrix;

    // set counters to 0
    memset( elements_assigned_to_process_counter, 0, n_processes*sizeof(int) );

    int i, p;
    for (i=0; i<matrix.n_elements; i++){
        int found = 0;
        for (p=0; p<n_processes && !found; p++){
            // if process P must handle element I
            if (check_if_inside_submatrix(partition->assignments[p], matrix.rows[i], matrix.cols[i])){
                found=1;
                // increase the counter of the elements handled by P
                elements_assigned_to_process_counter[p]++;
            }
        }
    }

    // MPI datatype
    build_mpi_tuple(&mpi_tuple);

    // send infos to other processes
    for (p=0; p<n_processes; p++){
        if (partition->n_elements_biggest_partition < elements_assigned_to_process_counter[p])
            partition->n_elements_biggest_partition = elements_assigned_to_process_counter[p];
        
        // send some informations to the specific process (so that he can allocate the correct space in memory)
        int infos[7];
        infos[0] = partition->assignments[p].start_row;
        infos[1] = partition->assignments[p].stop_row;
        infos[2] = partition->assignments[p].start_col;
        infos[3] = partition->assignments[p].stop_col;
        infos[4] = elements_assigned_to_process_counter[p];
        infos[5] = partition->assignments[p].col_responsible;
        infos[6] = partition->assignments[p].row_responsible;
  
        // create a countiguous block of memory to send the data
        struct mpi_matrix_element* mpi_data = malloc(sizeof(struct mpi_matrix_element)*elements_assigned_to_process_counter[p]);
        int count = 0;
        for (i=0; i<matrix.n_elements; i++){
            if (check_if_inside_submatrix(partition->assignments[p], matrix.rows[i], matrix.cols[i])){
                mpi_data[count].row = matrix.rows[i];
                mpi_data[count].col = matrix.cols[i];
                mpi_data[count].val = matrix.values[i];
                count++;
            }
        }
        if (p==0){
            process_0_submatrix = create_empty_submatrix(infos); //// TO BE DONE
            process_0_submatrix.elements = mpi_data;
        }else{
            MPI_Send( infos , 7 , MPI_INT , p , 0 , MPI_COMM_WORLD);
            MPI_Send( mpi_data , elements_assigned_to_process_counter[p] , mpi_tuple , p , 0 , MPI_COMM_WORLD);
            free(mpi_data);
        }
    }    
    MPI_Type_free(mpi_tuple);
    return process_0_submatrix;
}

submatrix create_empty_submatrix(int* infos){
    submatrix my_submatrix;

    // general infos
    my_submatrix.start_row = infos[0];
    my_submatrix.stop_row = infos[1];
    my_submatrix.start_col = infos[2];
    my_submatrix.stop_col = infos[3];
    my_submatrix.n_elements = infos[4];
    my_submatrix.col_responsible = infos[5];
    my_submatrix.row_responsible = infos[6];
    
    // elements of array
    my_submatrix.elements = malloc(sizeof(struct mpi_matrix_element)*infos[4]);

    return my_submatrix;
}

void clean_submatrix(submatrix *submatrix) {
    free(submatrix->elements);
    submatrix->elements = NULL;
}

submatrix wait_for_sparse_matrix(){
    int infos[7];
    MPI_Status status;
    submatrix my_submatrix;
    MPI_Datatype mpi_tuple;

    // initialize submatrix
    MPI_Recv( infos , 7 , MPI_INT , 0 , 0 , MPI_COMM_WORLD , &status);
    my_submatrix = create_empty_submatrix(infos);

    // add elements to the queue
    build_mpi_tuple(&mpi_tuple);
    MPI_Recv( my_submatrix.elements , my_submatrix.n_elements , mpi_tuple , 0 , 0 , MPI_COMM_WORLD , &status);

    int i;
    for (i = 0; i < my_submatrix.n_elements; i++) {
        my_submatrix.elements[i].row -= my_submatrix.start_row;
        my_submatrix.elements[i].col -= my_submatrix.start_col;
    }
    MPI_Type_free(mpi_tuple);
    return my_submatrix;
}


double_dense_matrix sum_submatrix_along_cols(const submatrix submatrix) {
    double_dense_matrix result_matrix = create_double_dense_matrix(1, submatrix.stop_col - submatrix.start_col);
    int i;
    for (i = 0; i < submatrix.stop_col - submatrix.start_col; i++) {
        result_matrix.matrix[i] = 0.0;
    }
    for (i = 0; i < submatrix.n_elements; i++) {
        result_matrix.matrix[submatrix.elements[i].col] += submatrix.elements[i].val;
    }
    return result_matrix;
}

double_dense_matrix sum_submatrix_along_rows(const submatrix submatrix) {
    double_dense_matrix result_matrix = create_double_dense_matrix(submatrix.stop_row - submatrix.start_row, 1);
    int i;
    for (i = 0; i < submatrix.stop_row - submatrix.start_row; i++) {
        result_matrix.matrix[i] = 0.0;
    }
    for (i = 0; i < submatrix.n_elements; i++) {
        result_matrix.matrix[submatrix.elements[i].row] += submatrix.elements[i].val;
    }
    return result_matrix;
}

void multiply_coefficient_by_cols(const double_dense_matrix alfa_i, submatrix working_submatrix) {
    assert(alfa_i.n_rows == 1);
    int i;
    for (i = 0; i < working_submatrix.n_elements; i++) {
        working_submatrix.elements[i].val *= alfa_i.matrix[working_submatrix.elements[i].col];
    }
}

void multiply_coefficient_by_rows(const double_dense_matrix alfa_i, submatrix working_submatrix) {
    assert(alfa_i.n_cols == 1);
    int i;
    for (i = 0; i < working_submatrix.n_elements; i++) {
        working_submatrix.elements[i].val *= alfa_i.matrix[working_submatrix.elements[i].row];
    }
}

submatrix clone_submatrix(const submatrix original_submatrix){
    submatrix clone;
    memcpy(&clone, &original_submatrix, sizeof(submatrix));
    clone.elements = malloc(sizeof(struct mpi_matrix_element)*(original_submatrix.n_elements+1));
    memcpy(clone.elements, original_submatrix.elements, sizeof(struct mpi_matrix_element)*original_submatrix.n_elements);
    clone.elements[clone.n_elements].row = -1; /// this one is needed when we send back the matrix --> know the length;
    return clone;
}

// sending submatrices back to process 0
void send_submatrices(submatrix working_submatrix){
    MPI_Datatype mpi_tuple;
    build_mpi_tuple(&mpi_tuple);
    int i;
    for (i = 0; i < working_submatrix.n_elements; i++) {
        working_submatrix.elements[i].row += working_submatrix.start_row;
        working_submatrix.elements[i].col += working_submatrix.start_col;
    }
    MPI_Send( working_submatrix.elements , working_submatrix.n_elements+1 , mpi_tuple , 0 , 0 , MPI_COMM_WORLD);
    MPI_Type_free(mpi_tuple);
}

double_sparse_matrix group_submatrices(submatrix working_submatrix, double_sparse_matrix* original_matrix, const int n_processes, const int n_elements_biggest_sumbatrix){
    double_sparse_matrix results = create_double_sparse_matrix(original_matrix->n_rows, original_matrix->n_cols, original_matrix->n_elements);
    struct mpi_matrix_element* buffer = malloc((n_elements_biggest_sumbatrix+1) * sizeof(struct mpi_matrix_element));
    int index = 0;
    MPI_Status status;
    MPI_Datatype mpi_tuple;
    build_mpi_tuple(&mpi_tuple);

    int j;
    // add the results of the local working_submatrix
    for (j=0; j<working_submatrix.n_elements; j++){
        results.rows[index] = working_submatrix.elements[j].row;
        results.cols[index] = working_submatrix.elements[j].col;
        results.values[index] = working_submatrix.elements[j].val;
        index++;
    }

    int i;
    // receive results from the other processes
    for (i = 1; i< n_processes; i++){
        j = 0;
        MPI_Recv( buffer , n_elements_biggest_sumbatrix+1 , mpi_tuple , i , 0 , MPI_COMM_WORLD , &status);

        while (buffer[j].row!=-1){
            results.rows[index] = buffer[j].row;
            results.cols[index] = buffer[j].col;
            results.values[index] = buffer[j].val;
            j++;
            index++;
        }
        // add element to results;
    }
    free(buffer);
    MPI_Type_free(mpi_tuple);
    return results;
}


void util_print_submatrix(const submatrix submatrix){
    printf("n_elements=%d   starts=[%d,%d] stop=[%d,%d]\n", submatrix.n_elements, submatrix.start_row, submatrix.start_col, submatrix.stop_row, submatrix.stop_col);
    int i;
    for (i=0; i<submatrix.n_elements; i++){
        printf("- (%d,%d): %lf\n", submatrix.elements[i].row, submatrix.elements[i].col, submatrix.elements[i].val);
    }
    printf("\n");
}