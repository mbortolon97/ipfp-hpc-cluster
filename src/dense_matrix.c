#include "dense_matrix.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>


/**
 * This function return a new dense matrix with the given size
 **/ 
double_dense_matrix create_double_dense_matrix(int n_rows, int n_cols) {
    double_dense_matrix matrix;

    matrix.n_rows = n_rows;
    matrix.n_cols = n_cols;
    matrix.matrix = malloc(n_rows * n_cols * sizeof(double));

    return matrix;
}

/**
 * This function free a dense matrix
**/
void clean_double_dense_matrix(double_dense_matrix* matrix) {
    free(matrix->matrix);
    matrix->matrix = NULL;
}

double_dense_matrix get_row_from_dense(double_dense_matrix matrix, int row) {
    int i;
    double_dense_matrix result_matrix = create_double_dense_matrix(1, matrix.n_cols);
    for (i = 0; i < matrix.n_cols; i++) {
        result_matrix.matrix[i] = matrix.matrix[row * matrix.n_cols + i];
    }
    return result_matrix;
}

double_dense_matrix get_col_from_dense(double_dense_matrix matrix, int col) {
    int i;
    double_dense_matrix result_matrix = create_double_dense_matrix(matrix.n_rows, 1);
    for (i = 0; i < matrix.n_rows; i++) {
        result_matrix.matrix[i] = matrix.matrix[i * matrix.n_cols + col];
    }
    return result_matrix;
}

void set_to_one_less_than_epsilon(double_dense_matrix matrix) {
    int i, j;
    for (i = 0; i < matrix.n_rows; i++) {
        for (j = 0; j < matrix.n_cols; j++) {
            if (matrix.matrix[i * matrix.n_cols + j] < FLT_EPSILON) {
                matrix.matrix[i * matrix.n_cols + j] = 1.0;
            }
            // matrix.matrix[i * matrix.n_cols + j] = (matrix.matrix[i * matrix.n_cols + j] < DBL_EPSILON) * 1.0 + (matrix.matrix[i * matrix.n_cols + j] > DBL_EPSILON) * matrix.matrix[i * matrix.n_cols + j];
        }
    }
}

/**
 * This function divide two dense array element wise
**/
double_dense_matrix elementwise_division(const double_dense_matrix matrix1, const double_dense_matrix matrix2) {
    assert(matrix1.n_rows == matrix2.n_rows);
    assert(matrix1.n_cols == matrix2.n_cols);

    double_dense_matrix result_matrix = create_double_dense_matrix(matrix1.n_rows, matrix1.n_cols);

    int i, j;
    for (i = 0; i < result_matrix.n_rows; i++) {
        for (j = 0; j < result_matrix.n_cols; j++) {
            result_matrix.matrix[i * result_matrix.n_cols + j] = matrix1.matrix[i * matrix1.n_cols + j] / matrix2.matrix[i * matrix2.n_cols + j];
        }
    }

    return result_matrix;
}

void print_dense_matrix(const double_dense_matrix matrix) {
    int i, j;
    for (i = 0; i < matrix.n_rows; i++) {
        for (j = 0; j < matrix.n_cols; j++) {
            printf("%e ", matrix.matrix[i * matrix.n_cols + j]);
        }
        printf("\n");
    }
}

/**
 * This function load from a file a double dense matrix
**/
double_dense_matrix load_double_dense_matrix(const char *filename) {
    int n_rows, n_cols, i;

    // load file
    FILE *in_file = fopen(filename, "r");

    // check file was opened correctly
    if (in_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    fscanf(in_file, "%d %d", &n_rows, &n_cols);

    double_dense_matrix loaded_matrix = create_double_dense_matrix(n_rows, n_cols);

    // read each element of the matrix
    for (i = 0; i < n_rows * n_cols; i++)
    {
        int row;
        int col;
        double value;
        fscanf(in_file, "%d %d %lf", &row, &col, &value);

        loaded_matrix.matrix[row * n_cols + col] = value;
    }
    fclose(in_file);

    return loaded_matrix;
}