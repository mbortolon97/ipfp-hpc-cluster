//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "input.h"

int *get_swapped_array(int len)
{
    int *array, tmp, j;
    array = (int *)malloc(len * sizeof(int));

    // initialize array
    for (int i = 0; i < len; i++)
        array[i] = i;

    // shuffle it
    for (int i = 0; i < len; i++)
    {
        j = rand() % len;
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
    return array;
}

struct sparse_matrix *load_sparse_matrix()
{
    char line[64], *ignore;
    int n_rows, n_cols, n_elements, i;
    struct sparse_matrix_element *elements;
    struct sparse_matrix *matrix;

    // load file
    FILE *in_file = fopen("matrix.txt", "r");

    // check file was opened correctly
    if (in_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    // get infos from first line of the matrix
    ignore = fgets(line, 64, in_file);
    n_rows = atoi(strtok(line, " "));
    n_cols = atoi(strtok(NULL, " "));
    n_elements = atoi(strtok(NULL, " "));

    // allocate space for the matrix
    elements = (struct sparse_matrix_element *)malloc(n_elements * sizeof(struct sparse_matrix_element));

    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        ignore = fgets(line, 64, in_file);
        elements[i].row = atof(strtok(line, " "));
        elements[i].col = atof(strtok(NULL, " "));
        elements[i].value = atof(strtok(NULL, " "));
    }
    fclose(in_file);

    // create matrix
    matrix = (struct sparse_matrix *)malloc(sizeof(struct sparse_matrix));
    matrix->n_rows = n_rows;
    matrix->n_cols = n_cols;
    matrix->n_elements = n_elements;
    matrix->elements = elements;

    // shuffle rows and columns
    matrix->row_permutations = NULL;
    matrix->col_permutations = NULL;

    return matrix;
}

void save_sparse_matrix(struct sparse_matrix *matrix)
{
    // open file
    char line[64];
    FILE *out_file = fopen("out_matrix.txt", "w+");

    // check file was opened correctly
    if (out_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1);
    }

    // write the first line
    sprintf(line, "%d %d %d\n", matrix->n_rows, matrix->n_cols, matrix->n_elements);
    fputs(line, out_file);

    // write elements
    for (int i = 0; i < matrix->n_elements; i++)
    {
        sprintf(line, "%d %d %f\n", matrix->elements[i].row, matrix->elements[i].col, matrix->elements[i].value);
        fputs(line, out_file);
    }
}

void permute_sparse_matrix(struct sparse_matrix *matrix)
{
    if (matrix->row_permutations != NULL)
        free(matrix->row_permutations);
    if (matrix->col_permutations != NULL)
        free(matrix->col_permutations);

    matrix->row_permutations = get_swapped_array(matrix->n_rows);
    matrix->col_permutations = get_swapped_array(matrix->n_cols);
}

void free_sparse_matrix(struct sparse_matrix *matrix)
{
    free(matrix->row_permutations);
    free(matrix->col_permutations);
    free(matrix->elements);
    free(matrix);
}