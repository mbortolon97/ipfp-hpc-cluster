#ifndef INPUT
#define INPUT

/**********************************
 * format of a matrix.txt file
 *
 * n_rows n_cols n_elements
 * x1 y1 element_value_1
 * x2 y2 element_value_2
 * x3 y3 element_value_3
 * .....
 * 
**********************************/

/* one element of the matrix */
struct sparse_matrix_element
{
    int row;
    int col;
    float value;
};

/* a struct representing a sparse matrix */
struct sparse_matrix
{
    int n_rows;
    int n_cols;
    int n_elements;
    int *row_permutations;
    int *col_permutations;
    struct sparse_matrix_element *elements;
};

/* load a file matrix.txt which contains a sparse matrix */
/* returns a pointer to a sparse_matrix */
struct sparse_matrix *load_sparse_matrix();

/* permute rows and columns */
void permute_sparse_matrix(struct sparse_matrix *matrix);

/* given a matrix, stores it in a out_matrix.txt file */
void save_sparse_matrix(struct sparse_matrix *matrix);

/* free memory in heap dedicated to the matrix */
void free_sparse_matrix(struct sparse_matrix *matrix);

#endif