#ifndef __PERMUTATION_H
#define __PERMUTATION_H 1

typedef struct sparse_matrix_permutation_struct
{
    int n_rows;
    int n_cols;
    int[] row_permutations;
    int[] col_permutations;
    int[] inverse_row_permutations;
    int[] inverse_col_permutations;
} sparse_matrix_permutation;

/**
 * This function return a new sparse matrix permutation for the given matrix
 **/
sparse_matrix_permutation create_sparse_matrix_random_permutation(const sparse_matrix matrix);

/**
 * This function free a sparse matrix permutation
**/
void clean_sparse_matrix_permutation(sparse_matrix_permutation* permutation);

/**
 * Given a sparse matrix it create a new sparse matrix permutated based on the given map
 **/
sparse_matrix permutate_double_sparse_matrix(const sparse_matrix_permutation permutation, const sparse_matrix matrix);

/**
 * Permutate a given dense matrix based on column of the sparse permutation
 **/
dense_matrix  permutate_dense_matrix_along_columns(const sparse_matrix_permutation permutation, const dense_matrix matrix);

/**
 * Permutate a given dense matrix based on row of the sparse permutation
 **/
dense_matrix  permutate_dense_matrix_along_rows(const sparse_matrix_permutation permutation, const dense_matrix matrix);


#endif