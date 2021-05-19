//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "matrix.h"
#include "input.h"

struct sparse_matrix *load_sparse_matrix()
{
    char line[64], *pch;
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
    fgets(line, 64, in_file);
    n_rows = atoi(strtok(line, " "));
    n_cols = atoi(strtok(NULL, " "));
    n_elements = atoi(strtok(NULL, " "));

    // allocate space for the matrix
    elements = (struct sparse_matrix_element *)malloc(n_elements * sizeof(struct sparse_matrix_element));

    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        fgets(line, 64, in_file);
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

/*
struct dense_matrix *load_dense_matrix()
{
    char line[64], *pch;
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
    fgets(line, 64, in_file);
    n_rows = atoi(strtok(line, " "));
    n_cols = atoi(strtok(NULL, " "));
    n_elements = atoi(strtok(NULL, " "));

    // allocate space for the matrix
    elements = (struct sparse_matrix_element *)malloc(n_elements * sizeof(struct sparse_matrix_element));

    // read each element of the matrix
    for (i = 0; i < n_elements; i++)
    {
        fgets(line, 64, in_file);
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

    return matrix;
}*/