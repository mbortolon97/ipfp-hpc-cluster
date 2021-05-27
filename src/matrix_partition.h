typedef struct submatrix_partition_struct {
    int subp_rows;
    int subp_cols;

    int[] process_row;
    int[] process_col;
    int[] rows_responsible;
    int[] cols_responsible;
} submatrix_partition;

submatrix_partition create_submatrix_partition(int n_processes, int n_rows, int n_cols);