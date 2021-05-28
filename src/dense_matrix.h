#ifndef __DENSE_MATRIX_H
#define __DENSE_MATRIX_H 1

typedef struct dense_matrix_struct
{
    int n_rows;
    int n_cols;
    double[] matrix;
} double_dense_matrix;

/**
 * This function return a new dense matrix with the given size
 **/ 
double_dense_matrix create_double_dense_matrix(int n_rows, int n_cols);

/**
 * This function free a dense matrix
**/
void clean_double_dense_matrix(double_dense_matrix* matrix);

/**
 * This function load from a file a double dense matrix
**/
double_dense_matrix load_double_dense_matrix(const char *filename);

/**
 * This function extract as a dense matrix a column
 **/ 
double_dense_matrix get_col_as_dense(const double_sparse_matrix matrix, int col);

#endif