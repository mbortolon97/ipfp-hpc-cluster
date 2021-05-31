#ifndef __DENSE_MATRIX_H
#define __DENSE_MATRIX_H 1

typedef struct dense_matrix_struct
{
    int n_rows;
    int n_cols;
    double* matrix;
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
 * This function get a row as a dense matrix
 **/
double_dense_matrix get_row_from_dense(double_dense_matrix matrix, int row);

/**
 * This function set to one all the value below epsilon, work in place
 **/
void set_to_one_less_than_epsilon(double_dense_matrix matrix);

/**
 * This function divide two dense array element wise
**/
double_dense_matrix elementwise_division(const double_dense_matrix matrix1, const double_dense_matrix matrix2);

void print_dense_matrix(const double_dense_matrix matrix);

#endif