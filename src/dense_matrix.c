#include "dense_matrix.h"

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