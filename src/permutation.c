#include "permutation.h"
#include "sparse_matrix.h"

typedef struct sparse_matrix_permutation_struct
{
    int n_rows;
    int n_cols;
    int[] row_permutations;
    int[] col_permutations;
    int[] inverse_row_permutations;
    int[] inverse_col_permutations;
} sparse_matrix_permutation;

sparse_matrix_permutation create_sparse_matrix_random_permutation(const sparse_matrix matrix) {
    sparse_matrix_permutation permutation;
    permutation.n_rows = matrix.n_rows;
    permutation.n_cols = matrix.n_cols;
    permutation.row_permutations = malloc(matrix.n_rows * sizeof(int));
    permutation.col_permutations = malloc(matrix.n_cols * sizeof(int));
    permutation.inverse_row_permutations = malloc(matrix.n_rows * sizeof(int));
    permutation.inverse_col_permutations = malloc(matrix.n_cols * sizeof(int));

    int i, j;
    for (i = 0; i < matrix.n_rows; i++) {
        j = rand() % matrix.n_rows;
        permutation.row_permutations[i] = j;
        permutation.inverse_row_permutations[j] = i;
    }
    for (i = 0; i < matrix.n_cols; i++) {
        j = rand() % matrix.n_cols;
        permutation.col_permutations[i] = j;
        permutation.inverse_col_permutations[j] = i;
    }
}

void clean_sparse_matrix_permutation(sparse_matrix_permutation* permutation) {
    free(permutation->row_permutations);
    free(permutation->col_permutations);
    free(permutation->inverse_row_permutations);
    free(permutation->inverse_col_permutations);
    permutation->row_permutations = NULL;
    permutation->col_permutations = NULL;
    permutation->inverse_row_permutations = NULL;
    permutation->inverse_col_permutations = NULL;
}

sparse_matrix permutate_double_sparse_matrix(const sparse_matrix_permutation permutation, const sparse_matrix matrix) {
    sparse_matrix permutated_matrix = create_double_sparse_matrix(matrix.n_rows, matrix.n_cols, matrix.n_elements);
    int i;
    for (i = 0; i < matrix.n_elements; i++) {
        permutated_matrix.rows[i] = permutation.row_permutations[matrix.rows[i]];
        permutated_matrix.cols[i] = permutation.col_permutations[matrix.cols[i]];
        permutated_matrix.values[i] = matrix.values[i];
    }
    return permutated_matrix;
}

dense_matrix permutate_dense_matrix_along_columns(const sparse_matrix_permutation permutation, const dense_matrix matrix) {
    dense_matrix permutated_matrix = create_double_dense_matrix(matrix.n_rows, matrix.n_cols);

    int i, j;
    for (i = 0; i < matrix.n_rows; i++) {
        for (j = 0; j < matrix.n_cols; j++) {
            permutated_matrix.matrix[i * matrix.n_cols + j] = matrix.matrix[i * matrix.n_cols + permutation.col_permutations[j]];
        }
    }

    return permutated_matrix;
}

dense_matrix permutate_dense_matrix_along_rows(const sparse_matrix_permutation permutation, const dense_matrix matrix) {
    dense_matrix permutated_matrix = create_double_dense_matrix(matrix.n_rows, matrix.n_cols);

    int i, j;
    for (i = 0; i < matrix.n_rows; i++) {
        for (j = 0; j < matrix.n_cols; j++) {
            permutated_matrix.matrix[i * matrix.n_cols + j] = matrix.matrix[permutation.row_permutations[i] * matrix.n_cols + j];
        }
    }

    return permutated_matrix;
}