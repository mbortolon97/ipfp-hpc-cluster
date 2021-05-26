#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "input.h"

struct factor_queue_element {
    int factor;
    factor_queue_element* next;
};

// given a number to factorize (n), returns an array of factors (res) and its length (len)
struct factor_queue_element* factorize(int n, int *res, int *len)
{
    int factor = 2;
    *len = 0;

    while (n > 1)
    {
        while (n % factor != 0)
        {
            factor++;
        }
        res[*len] = factor;
        len[0] = len[0] + 1;
        n = n / factor;
    }
}

// uses heuristic to decide how to better split the matrix
void number_of_submatrices(int n_processes, int n_rows, int n_cols, int *subp_rows, int *subp_cols)
{
    int factors[32], len;
    float ideal_rateo, actual_rateo;

    // try two different factorizations
    factorize(n_processes, factors, &len);

    // if someone uses a prime number as number of processes....
    if (len == 1 && n_processes > 3)
        factorize(n_processes - 1, factors, &len);

    // calculating how to divide the matrix
    ideal_rateo = ((double)n_rows) / (double)n_cols;
    *subp_rows = 1;
    *subp_cols = 1;
    for (int i = len - 1; i >= 0; i--)
    {
        actual_rateo = ((double)*subp_rows) / (double)*subp_cols;
        if (actual_rateo >= ideal_rateo)
            // too rows --> add columns
            *subp_cols *= factors[i];
        else
            // too columns --> add rows
            *subp_rows *= factors[i];
    }
}