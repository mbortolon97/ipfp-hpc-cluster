#ifndef SMALLMAT
#define SMALLMAT

struct queue
{
    float value;
    struct queue *next;
};

struct queue_matrix
{
    int n_rows;
    int n_cols;
    struct queue **row_index;
    struct queue **col_index;
};

struct queue_matrix *make_matrix(int n_rows, int n_cols); /// to do --> take as input the stuff received from mpi

#endif