#include "dense_matrix.h"
#include <assert.h>

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
    double_dense_matrix result_matrix = create_double_dense_matrix(1, matrix.n_cols);
    for (i = 0; i < matrix.n_cols; i++) {
        result_matrix.matrix[i] = matrix.matrix[row * n_cols + i];
    }
    return result_matrix;
}

void set_to_one_less_than_epsilon(double_dense_matrix matrix) {
    int i, j;
    for (i = 0; i < matrix.n_rows; i++) {
        for (j = 0; j < matrix.n_cols; j++) {
            matrix.matrix[i * n_cols + j] = matrix.matrix[i * n_cols + j] < DBL_EPSILON * 1.0 + matrix.matrix[i * n_cols + j] > DBL_EPSILON * matrix.matrix[i * n_cols + j];
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
    for (i = 0; i < matrix.n_rows; i++) {
        for (j = 0; j < matrix.n_cols; j++) {
            result_matrix.matrix[i * n_cols + j] = matrix1.matrix[i * n_cols + j] / matrix2.matrix[i * n_cols + j];
        }
    }
    
    return result_matrix;
}

/**
 * This function load from a file a double dense matrix
**/
struct double_dense_matrix load_double_dense_matrix(const char *filename) {
    int n_rows, n_cols;

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

    double* pointer = (double*)matrix.matrix
    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        int row;
        int col;
        double value;
        fscanf(in_file, "%d %d %lf", &row, &col, &value);

        matrix.matrix[row * n_cols + col] = value;
    }
    fclose(in_file);

    return loaded_matrix;
}