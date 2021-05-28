#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H 1

typedef struct sparse_matrix_struct
{
    int n_rows;
    int n_cols;
    int n_elements;
    int[] rows;
    int[] cols;
    double[] values;
} double_sparse_matrix;

/**
 * This function return a new sparse matrix with the given size
 **/ 
double_sparse_matrix create_double_sparse_matrix(int n_rows, int n_cols, int n_elements);

/**
 * This function free a sparse matrix
**/
void clean_sparse_matrix(double_sparse_matrix* matrix);

/**
 * This function load from a file a double sparse matrix
**/
struct double_sparse_matrix load_double_sparse_matrix(const char *filename);

/**
 * This function save to a file a double sparse matrix
**/
void save_double_sparse_matrix(const char *filename, double_sparse_matrix matrix);

#endif