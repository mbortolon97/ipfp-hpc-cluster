#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "input.h"

// given a number to factorize (n), returns an array of factors (res) and its length (len)
void factorize(int n, int *res, int *len)
{
    int factor = 2;
    *len = 0;

    while (n > 1)
    {
        while (n % factor != 0)
        {
            factor++;
        }
        res[*len] = factor;
        len[0] = len[0] + 1;
        n = n / factor;
    }
}

// uses heuristic to decide how to better split the matrix
void number_of_submatrices(int n_processes, int n_rows, int n_cols, int *subp_rows, int *subp_cols)
{
    int factors[32], len;
    float ideal_rateo, actual_rateo;

    // try two different factorizations
    factorize(n_processes, factors, &len);

    // if someone uses a prime number as number of processes....
    if (len == 1 && n_processes > 3)
        factorize(n_processes - 1, factors, &len);

    // calculating how to divide the matrix
    ideal_rateo = ((double)n_rows) / (double)n_cols;
    *subp_rows = 1;
    *subp_cols = 1;
    for (int i = len - 1; i >= 0; i--)
    {
        actual_rateo = ((double)*subp_rows) / (double)*subp_cols;
        if (actual_rateo >= ideal_rateo)
            // too rows --> add columns
            *subp_cols *= factors[i];
        else
            // too columns --> add rows
            *subp_rows *= factors[i];
    }
}

void build_mpi_tuple(MPI_Datatype *mpi_tuple)
{
    struct sparse_matrix_element element;
    // array of structure member sizes
    int blocklengths[3] = {1, 1, 1};

    // structure member types
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_FLOAT};

    // offset of structure members
    MPI_Aint offsets[3];
    MPI_Aint points[3];
    MPI_Get_address(&element.row, points);
    MPI_Get_address(&element.col, points + 1);
    MPI_Get_address(&element.value, points + 2);
    offsets[0] = 0;
    offsets[1] = points[1] - points[0];
    offsets[2] = points[2] - points[0];

    // create mpi struct
    MPI_Type_create_struct(3, blocklengths, offsets, types, mpi_tuple);
    MPI_Type_commit(mpi_tuple);
}

struct sparse_matrix_element *split_in_submatrices(struct sparse_matrix *matrix, int n_processes, int *elements_per_process)
{
    int i, j, responsible_process, offset, elements_partial_sum[n_processes];
    int subp_rows, subp_cols, rows_per_subprocess, cols_per_subprocess, max_;
    struct sparse_matrix_element element;
    struct sparse_matrix_element *sorted_matrix = malloc(sizeof(struct sparse_matrix_element) * matrix->n_elements);

    srand(time(NULL));

    // get how many  &subp_rows, &subp_cols
    number_of_submatrices(n_processes, matrix->n_rows, matrix->n_cols, &subp_rows, &subp_cols);
    rows_per_subprocess = (matrix->n_rows) / subp_rows;
    cols_per_subprocess = (matrix->n_cols) / subp_cols;

    // permute columns and compute how many  elements would follow in each submatrix
    do
    {
        permute_sparse_matrix(matrix);
        // initialize counters to 0
        for (int i = 0; i < n_processes; i++)
            elements_per_process[i] = 0;

        // find how many elements each process should handle
        for (int k = 0; k < matrix->n_elements; k++)
        {
            element = matrix->elements[k];
            i = matrix->row_permutations[element.row]; // premuted row
            j = matrix->col_permutations[element.col]; // premuted col
            responsible_process = (i / rows_per_subprocess) * subp_cols;
            responsible_process += (j / cols_per_subprocess);
            elements_per_process[responsible_process]++;
        }

        // compute worst case
        max_ = 0;
        for (int i = 0; i < n_processes; i++)
            if (max_ < elements_per_process[i])
                max_ = elements_per_process[i];

        // if worst case is too bad: retry
    } while (max_ / 2 > matrix->n_elements / n_processes);

    elements_partial_sum[0] = elements_per_process[0];
    for (int i = 1; i < n_processes; i++)
        elements_partial_sum[i] = elements_partial_sum[i - 1] + elements_per_process[i];

    for (int k = 0; k < matrix->n_elements; k++)
    {
        element = matrix->elements[k];
        i = matrix->row_permutations[element.row];
        j = matrix->col_permutations[element.col];
        responsible_process = (i / rows_per_subprocess) * subp_cols;
        responsible_process += (j / cols_per_subprocess);
        offset = --elements_partial_sum[responsible_process];
        sorted_matrix[offset].row = i;
        sorted_matrix[offset].col = j;
        sorted_matrix[offset].value = element.value;
    }

    return sorted_matrix;
}