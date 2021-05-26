#include "sparse_matrix.h"
#include <stdio.h>

sparse_matrix create_double_sparse_matrix(int n_rows, int n_cols, int n_elements) {
    sparse_matrix matrix;
    
    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;
    matrix.n_elements = n_elements;
    matrix.rows = malloc(n_elements * sizeof(int));
    matrix.cols = malloc(n_elements * sizeof(int));
    matrix.values = malloc(n_elements * sizeof(double));
}

void clean_sparse_matrix(sparse_matrix* matrix) {
    free(matrix->rows);
    free(matrix->cols);
    free(matrix->values);
    matrix->rows = NULL;
    matrix->cols = NULL;
    matrix->values = NULL;
}

struct sparse_matrix* load_double_sparse_matrix(const char *filename)
{
    int n_rows, n_cols, n_elements, i;

    // load file
    FILE *in_file = fopen(filename, "r");

    // check file was opened correctly
    if (in_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    fscanf(in_file, "%d %d %d", &n_rows, &n_cols, &n_elements);

    loaded_matrix = create_sparse_matrix(n_rows, n_cols, n_elements);

    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        int row;
        int col;
        double value;
        fscanf(in_file, "%d %d %lf", &row, &col, &value);

        matrix.rows[i] = row;
        matrix.cols[i] = col;
        matrix.values[i] = value;
    }
    fclose(in_file);

    return loaded_matrix;
}

void save_double_sparse_matrix(const char *filename, sparse_matrix matrix)
{
    FILE *out_file = fopen(filename, "w+");

    // check file was opened correctly
    if (out_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    // write the first line
    fprintf(out_file, "%d %d %d\n", matrix.n_rows, matrix.n_cols, matrix.n_elements);

    // write elements
    for (int i = 0; i < matrix.n_elements; i++)
    {
        fprintf(out_file, "%d %d %lf\n", matrix.rows[i], matrix.cols[i], matrix.values[i]);
    }
}