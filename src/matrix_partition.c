#include "matrix_partition.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "input.h"

struct factor_list_element {
    int factor;
    factor_list_element* next;
    factor_list_element* previous;
};

struct factor_list {
    struct factor_list_element* last_element;
    struct factor_list_element* first_element;
}

// given a number to factorize (n), returns an array of factors (res) and its length (len)
struct factor_list factorize(int n) {
    int factor = 2;
    struct factor_list list;
    while (n > 1)
    {
        while (n % factor != 0)
        {
            factor++;
        }
        new_element = malloc(sizeof(struct factor_list_element));
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
    list->start_element = NULL;
    list->first_element = NULL;
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
    
    struct factor_list_element* current_element = list->last_element;
    while (current_element != NULL) {
        actual_rateo = ((double)*matrix_partition.subp_rows) / ((double)*matrix_partition.subp_cols);
        if (actual_rateo >= ideal_rateo) {
            // too rows --> add columns
            *(matrix_partition.subp_cols) *= current_element->factor;
        } else {
            // too columns --> add rows
            *(matrix_partition.subp_rows) *= current_element->factor;
        }
        
        current_element = current_element->previous;
    }

    clean_factor_list(&list)

    return matrix_partition;
}

submatrix_partition create_submatrix_partition(int n_processes, int n_rows, int n_cols) {
    submatrix_partition partition = number_of_submatrices(n_processes, n_rows, n_cols);

    int process_id = 0;
    for (int i = 0; i < partition.subp_rows; i++) {
        for (int j = 0; j < partition.subp_cols; j++) {
            process_id = partition
            process_id++;
        }
    }

    return partition;
}