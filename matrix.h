#ifndef MATRIX_H
#define MATRIX_H

#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

double count_residual (int n, int m, double *a, double *b, double *c, int k, int p, MPI_Comm comm);
void mpi_multiply_matrix_on_vector (double *a, double *b, double *res, int n, int m, int k, int p, MPI_Comm comm);
int read_matrix_col(double *a, int n, int m, int p, int k, double *buf, char *name, MPI_Comm comm);
void init_matrix(double *a, int s, int n, int m, int p, int k);
void init_vector(double *b, double *a, int n, int m, int p, int k, MPI_Comm comm);
void print_matrix_col(double *a, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm comm);
void copy_matrix(double *x, double *a, int n);
void copy_vector(double *c, double *b, int n);
double f(int s, int n, int i, int j);
double norma_vector(double *a, int n);
double norma_matrix(double *a, int m, int n);
void multiply_matrix_and_vector(double *a, double *b, double *c, int m, int n);
void multiply_blocks_m(double *a, double *b, double *c, int m);
void subtract_vectors(double *b, double *c, int n);
void print_vector(double *b, int n, int r);
int get_reversed_matrix_gauss(int n, double *a, double *x, int *r, double eps);
int solve_jordan (double *a, double *a_buf, double *b, double norm, int n, int m, int k, int p, MPI_Comm comm);
double mpi_norma_matrix (double *a, int m, int n, int col_cnt, int k, int p, MPI_Comm comm);


enum RETURN_CODES
{
    SUCCESS = 0,
    ERROR_READ = -1,
    ERROR_OPEN = -2,
    SING_MATRIX = -3,
    ERROR_MEMORY = -4,
    ERROR = -5
};

#endif
