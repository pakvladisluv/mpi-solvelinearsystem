
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define MAX 4
#define FIVE 5
#define DECIMAL 10
#define GIGA_MODIFIER 1e9
#define NANO_MODIFIER 1e-9

double func(int k, int n, int i, int j);
int matrix_input(int argc, int n, int k, int r, int s, int tm, int z, double *qb, double *a,
            double *b, char *argv[], MPI_Comm newcom);
struct timespec diff(struct timespec start, struct timespec end);
int Cholesky(double *a, double *ar, double *q, double *x, double *d, double *b, int n, int r, int s,
       int tm, MPI_Comm newcom);
void residual(const double *a, const double *q, double *b, double *x, int n, int r,
         int s, int tm, MPI_Comm newcom);
void out(double *a, double *b, double *qb, int n, int m, int r, int s, int tm,
         MPI_Comm newcom);
int input(char *t);
