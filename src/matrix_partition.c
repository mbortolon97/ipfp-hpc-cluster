#include "matrix_partition.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct factor_list_element {
    int factor;
    struct factor_list_element* next;
    struct factor_list_element* previous;
};

struct factor_list {
    struct factor_list_element* last_element;
    struct factor_list_element* first_element;
};


// returns true for the last process if it has to return (not used)
bool check_number_of_processes(int* world_size, int* world_rank){
    if (world_size[0]<5) return false;
    int half = (*world_size) / 2;

    int i;
    // if is not a prime number return 0
    for (i=2; i<half; i++)
        if (world_size[0]%i == 0) return false;
    
    // else decrease world size
    *world_size = (*world_size)-1;
    return *world_size == *world_rank;
}


// given a number to factorize (n), returns an array of factors (res) and its length (len)
struct factor_list factorize(int n) {
    int factor = 2;
    struct factor_list list;
    list.last_element = NULL;
    list.first_element = NULL;
    if (n == 1) {
        struct factor_list_element* new_element = malloc(sizeof(struct factor_list_element));
        new_element->next = NULL;
        new_element->factor = 1;
        new_element->previous = NULL;
        list.last_element = new_element;
        list.first_element = new_element;
        return list;
    }
    while (n > 1)
    {
        while (n % factor != 0)
        {
            factor++;
        }
        struct factor_list_element* new_element = malloc(sizeof(struct factor_list_element));
        new_element->next = NULL;
        new_element->factor = factor;
        new_element->previous = list.last_element;
        if (list.last_element == NULL) {
            list.first_element = new_element;
        } else {
            list.last_element->next = new_element;
        }
        list.last_element = new_element;
        n = n / factor;
    }
    return list;
}

void clean_factor_list(struct factor_list* list) {
    struct factor_list_element* next_element = list->first_element;
    while (next_element->next != NULL) {
        struct factor_list_element* tmp = next_element->next;
        free(next_element);
        next_element = tmp;
    }
    free(next_element);
    list->first_element = NULL;
    list->last_element = NULL;
}

// uses heuristic to decide how to better split the matrix
submatrix_partition number_of_submatrices(int n_processes, int n_rows, int n_cols)
{
    double ideal_rateo, actual_rateo;

    // try two different factorizations
    struct factor_list list = factorize(n_processes);

    // if someone uses a prime number as number of processes....
    if (list.first_element->next == NULL && n_processes > 3) {
        list = factorize(n_processes - 1);
    }
    

    // calculating how to divide the matrix
    ideal_rateo = ((double)n_rows) / (double)n_cols;
    submatrix_partition matrix_partition;
    matrix_partition.subp_rows = 1;
    matrix_partition.subp_cols = 1;
    
    struct factor_list_element* current_element = list.last_element;
    while (current_element != NULL) {
        actual_rateo = ((double)matrix_partition.subp_rows) / ((double)matrix_partition.subp_cols);
        if (actual_rateo >= ideal_rateo) {
            // too rows --> add columns
            matrix_partition.subp_cols *= current_element->factor;
        } else {
            // too columns --> add rows
            matrix_partition.subp_rows *= current_element->factor;
        }
        
        current_element = current_element->previous;
    }

    clean_factor_list(&list);

    return matrix_partition;
}

submatrix_partition create_submatrix_partition(int n_processes, int n_rows, int n_cols) {
    submatrix_partition partition = number_of_submatrices(n_processes, n_rows, n_cols);

    partition.assignments = malloc(n_processes * sizeof(submatrix_assignment));
    partition.col_master = malloc(partition.subp_cols * sizeof(submatrix_col_master));
    partition.row_master = malloc(partition.subp_rows * sizeof(submatrix_row_master));
    partition.n_elements_biggest_partition = 0;

    int process_id = 0;
    int n_rows_per_process = n_rows / partition.subp_rows;
    int n_cols_per_process = n_cols / partition.subp_cols;
    
    int i, j;
    for (i = 0; i < partition.subp_rows; i++) {
        partition.row_master[i].start_row = n_rows_per_process * i;
        partition.row_master[i].stop_row = n_rows_per_process * (i + 1);
        if (i == partition.subp_rows - 1) {
            partition.row_master[i].stop_row += n_rows - (partition.subp_rows * n_rows_per_process);
        }
        partition.row_master[i].row_master_process_id = -1;
        partition.row_master[i].num_subprocesses = 0;
        partition.row_master[i].row_processes_id = malloc((partition.subp_cols - 1) * sizeof(int));
    }
    for (j = 0; j < partition.subp_cols; j++) {
        partition.col_master[j].start_col = n_cols_per_process * j;
        partition.col_master[j].stop_col = n_cols_per_process * (j + 1);
        if (j == partition.subp_cols - 1) {
            partition.col_master[j].stop_col += n_cols - (partition.subp_cols * n_cols_per_process);
        }
        partition.col_master[j].col_master_process_id = -1;
        partition.col_master[j].num_subprocesses = 0;
        partition.col_master[j].col_processes_id = malloc((partition.subp_rows - 1) * sizeof(int));
    }
    for (i = 0; i < partition.subp_rows; i++) {

        for (j = 0; j < partition.subp_cols; j++) {
            if (partition.row_master[i].row_master_process_id == -1) {
                partition.row_master[i].row_master_process_id = process_id;
            } else {
                partition.row_master[i].row_processes_id[partition.row_master[i].num_subprocesses] = process_id;
                partition.row_master[i].num_subprocesses++;
            }
            if (partition.col_master[j].col_master_process_id == -1) {
                partition.col_master[j].col_master_process_id = process_id;
            } else {
                partition.col_master[j].col_processes_id[partition.col_master[j].num_subprocesses] = process_id;
                partition.col_master[j].num_subprocesses++;
            }
            partition.assignments[process_id].start_row = n_rows_per_process * i;
            partition.assignments[process_id].stop_row = n_rows_per_process * (i + 1);
            if (i == partition.subp_rows - 1) {
                partition.assignments[process_id].stop_row += n_rows - (partition.subp_rows * n_rows_per_process);
            }
            partition.assignments[process_id].start_col = n_cols_per_process * j;
            partition.assignments[process_id].stop_col = n_cols_per_process * (j + 1);
            if (j == partition.subp_cols - 1) {
                partition.assignments[process_id].stop_col += n_cols - (partition.subp_cols * n_cols_per_process);
            }

            partition.assignments[process_id].col_responsible = partition.col_master[j].col_master_process_id;
            partition.assignments[process_id].row_responsible = partition.row_master[i].row_master_process_id;
            process_id++;
        }
    }

    return partition;
}

void print_submatrix(const submatrix_partition partition) {
    printf("partition.subp_rows %d\n", partition.subp_rows);
    printf("partition.subp_cols %d\n", partition.subp_cols);
    for (i = 0; i < partition.subp_rows; i++) {
        printf("row: %d partition.row_master[i].start_row %d\n", i, partition.row_master[i].start_row);
        printf("row: %d partition.row_master[i].stop_row %d\n", i, partition.row_master[i].stop_row);
        printf("row: %d partition.row_master[i].num_subprocesses %d\n", i, partition.row_master[i].num_subprocesses);
        printf("row: %d [", i);
        for (int i = 0; i < partition.row_master[i].num_subprocesses; i++) {
            printf("%d, ", partition.row_master[i].row_processes_id[i]);
        }
        printf("]\n");
        printf("row: %d partition.row_master[i].row_master_process_id %d\n", i, partition.row_master[i].row_master_process_id);
    }
    for (i = 0; i < partition.subp_cols; i++) {
        printf("row: %d partition.col_master[i].start_col %d\n", i, partition.col_master[i].start_col);
        printf("row: %d partition.col_master[i].stop_col %d\n", i, partition.col_master[i].stop_col);
        printf("row: %d partition.col_master[i].num_subprocesses %d\n", i, partition.col_master[i].num_subprocesses);
        printf("row: %d [", i);
        for (int i = 0; i < partition.col_master[i].num_subprocesses; i++) {
            printf("%d, ", partition.col_master[i].col_processes_id[i]);
        }
        printf("]\n");
        printf("row: %d partition.col_master[i].col_master_process_id %d\n", i, partition.col_master[i].col_master_process_id);
    }
}

void clean_submatrix_partition(submatrix_partition* partition) {
    free(partition->assignments);
    int i;
    for (i = 0; i < partition->subp_rows; i++) {
        free(partition->row_master[i].row_processes_id);
    }
    for (i = 0; i < partition->subp_cols; i++) {
        free(partition->col_master[i].col_processes_id);
    }
    free(partition->col_master);
    free(partition->row_master);
}