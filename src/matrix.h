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
    sparse_matrix_element* elements;
}

struct dense_matrix
{
    int n_rows;
    int n_cols;
    float* matrix;
}

struct queue_matrix
{
    int n_rows;
    int n_cols;
    float* matrix;
}