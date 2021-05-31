#include "sparse_matrix.h"
#include <stdio.h>
#include <stdlib.h>

double_sparse_matrix create_double_sparse_matrix(int n_rows, int n_cols, int n_elements) {
    double_sparse_matrix matrix;
    
    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;
    matrix.n_elements = n_elements;
    matrix.rows = malloc(n_elements * sizeof(int));
    matrix.cols = malloc(n_elements * sizeof(int));
    matrix.values = malloc(n_elements * sizeof(double));

    return matrix;
}

void clean_sparse_matrix(double_sparse_matrix* matrix) {
    free(matrix->rows);
    free(matrix->cols);
    free(matrix->values);
    matrix->rows = NULL;
    matrix->cols = NULL;
    matrix->values = NULL;
}

double_dense_matrix get_col_as_dense(const double_sparse_matrix matrix, int col) {
    double_dense_matrix result_matrix = create_double_dense_matrix(matrix.n_rows, 1);
    int i;
    for (int i = 0; i < matrix.n_rows; i++) {
        result_matrix.matrix[i] = 0.0;
    }
    for (i = 0; i < matrix.n_elements; i++) {
        if (matrix.cols[i] == col) {
            result_matrix.matrix[matrix.rows[i]] += matrix.values[i];
        }
    }
    return result_matrix;
}

double_dense_matrix get_row_as_dense(const double_sparse_matrix matrix, int row) {
    double_dense_matrix result_matrix = create_double_dense_matrix(1, matrix.n_cols);
    int i;
    for (int i = 0; i < matrix.n_cols; i++) {
        result_matrix.matrix[i] = 0.0;
    }
    for (i = 0; i < matrix.n_elements; i++) {
        if (matrix.rows[i] == row) {
            result_matrix.matrix[matrix.cols[i]] += matrix.values[i];
        }
    }
    return result_matrix;
}

double_sparse_matrix load_double_sparse_matrix(const char *filename)
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

    double_sparse_matrix loaded_matrix = create_double_sparse_matrix(n_rows, n_cols, n_elements);

    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        int row;
        int col;
        double value;
        fscanf(in_file, "%d %d %lf", &row, &col, &value);

        loaded_matrix.rows[i] = row;
        loaded_matrix.cols[i] = col;
        loaded_matrix.values[i] = value;
    }
    fclose(in_file);

    return loaded_matrix;
}

void save_double_sparse_matrix(const char *filename, double_sparse_matrix matrix)
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
        fprintf(out_file, "%d %d %e\n", matrix.rows[i], matrix.cols[i], matrix.values[i]);
    }
}