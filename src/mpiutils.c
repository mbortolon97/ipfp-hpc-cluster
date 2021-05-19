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
int number_of_submatrices(int n_processes, int n_rows, int n_cols, int *subp_rows, int *subp_cols)
{
    int factors[32], len;
    float ideal_rateo, actual_rateo, dropped = 0;

    // try two different factorizations
    factorize(n_processes, factors, &len);

    // if someone uses a prime number as number of processes....
    if (len == 1 && n_processes > 3)
    {
        factorize(n_processes - 1, factors, &len);
        dropped = 1;
    }

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
    return dropped;
}

void split_in_submatrices(struct sparse_matrix *matrix, int n_processes)
{
    int element_per_process[n_processes], i, j, responsible_process;
    int subp_rows, subp_cols, dropped, rows_per_subprocess, cols_per_subprocess, max_;
    struct sparse_matrix_element element;
    srand(time(NULL));

    // get how many  &subp_rows, &subp_cols
    dropped = number_of_submatrices(n_processes, matrix->n_rows, matrix->n_cols, &subp_rows, &subp_cols);
    rows_per_subprocess = (matrix->n_rows) / subp_rows;
    cols_per_subprocess = (matrix->n_cols) / subp_cols;

    // permute columns and compute how many  elements would follow in each submatrix
    do
    {
        permute_sparse_matrix(matrix);
        // initialize counters to 0
        for (int i = 0; i < n_processes; i++)
            element_per_process[i] = 0;

        // find how many elements each process should handle
        for (int k = 0; k < matrix->n_elements; k++)
        {
            element = matrix->elements[k];
            i = matrix->row_permutations[element.row]; // premuted row
            j = matrix->col_permutations[element.col]; // premuted col
            responsible_process = (i / rows_per_subprocess) * subp_cols;
            responsible_process += (j / cols_per_subprocess);
            element_per_process[responsible_process]++;
        }

        // compute worst case
        max_ = 0;
        for (int i = 0; i < n_processes; i++)
            if (max_ < element_per_process[i])
                max_ = element_per_process[i];

        // if worst case is too bad: retry
    } while (max_ / 2 > matrix->n_elements / n_processes);
}