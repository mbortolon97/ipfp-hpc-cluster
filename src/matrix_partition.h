typedef struct submatrix_assignment_struct {
    int start_row;
    int stop_row;
    int start_col;
    int stop_col;

    int col_responsible;
    int row_responsible;
} submatrix_assignment;

typedef struct submatrix_row_master_struct {
    int start_row;
    int stop_row;
    int row_master_process_id;
    int num_subprocesses;
    int[] row_processes_id;
} submatrix_row_master;

typedef struct submatrix_col_master_struct {
    int start_col;
    int stop_col;
    int col_master_process_id;
    int num_subprocesses;
    int[] col_processes_id;
} submatrix_col_master;

typedef struct submatrix_partition_struct {
    int subp_rows;
    int subp_cols;

    submatrix_assignment[] assignments;
    submatrix_col_master[] col_master;
    submatrix_row_master[] row_master;

    int[] process_row;
    int[] process_col;
    int[] rows_responsible;
    int[] cols_responsible;
} submatrix_partition;

submatrix_partition create_submatrix_partition(int n_processes, int n_rows, int n_cols);

void clean_submatrix_partition(submatrix_partition* partition);

uint8_t inline check_if_inside_submatrix(submatrix_assignment assignment, int row, int col);