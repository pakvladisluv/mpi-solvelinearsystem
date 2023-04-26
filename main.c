#include "func.h"

int main(int argc, char **argv)
{
    int r;
    int s;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
    MPI_Comm_size(MPI_COMM_WORLD, &s);
    if (argc < 4 || argc > FIVE) {
        printf("Неверное количество аргументов\n");
        return -1;
    }
    int n = input(argv[1]);
    int m = input(argv[2]);
    int k = input(argv[3]);
    int z = n % s;
    if (n <= 0 || m <= 0 || m > n || k < 0 || k > 4) {
        printf("некорректные значения аргументов\n");
        return -1;
    }
    if (k == 0 && argc == 4) {
        printf("при k == 0 не указан файл\n");
        return -1;
    }
    int q = 1;
    int tm = MPI_UNDEFINED;
    if (z > r) {
        tm = 0;
        q = 0;
    }
    q*=1;
    MPI_Comm newcom;
    MPI_Comm_split(MPI_COMM_WORLD, tm, r, &newcom);
    double *b = NULL;
    double *qb = NULL;
    if (r == 0) {
        b = (double *)malloc(n * sizeof(double));
    }
    qb = (double *)malloc(n * sizeof(double));
    double *a = (double *)malloc(n * (n / s + 1) * sizeof(double));
    if (!a) {
        free(a);
        free(qb);
        if (r == 0) {
            free(b);
        }
        return -2;
    }
    if (matrix_input(argc, n, k, r, s, tm, z, qb, a, b, argv, newcom) == -3) {
        if (r == 0) {
            free(b);
        }
        free(qb);
        free(a);
        return -3;
    }
    struct timespec time_start;
    struct timespec time_end;
    double *x = (double *)malloc(n * sizeof(double));
    double *d = (double *)malloc(n * sizeof(double));
    out(a, b, qb, n, m, r, s, tm, newcom);
    double *ar = (double *)malloc(n * (n / s + 1) * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n / s + 1; j++) {
            ar[i * (n / s + 1) + j] = 0;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(qb, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (r == 0) {
        clock_gettime(CLOCK_MONOTONIC, &time_start);
    }
    if (Cholesky(a, ar, qb, x, d, b, n, r, s, q, newcom) == -4) {
        free(a);
        free(ar);
        if (r == 0) {
            free(b);
        }
        free(qb);
        free(x);
        free(d);
        return -4;
    }
    if (r == 0) {
        clock_gettime(CLOCK_MONOTONIC, &time_end);
        time_end = diff(time_start, time_end);
        printf("Time: %lf s\n",
               (double)(time_end.tv_sec + time_end.tv_nsec * NANO_MODIFIER));
        printf("Solution: ");
        for (int i = 0; i < m; i++) {
            printf("%10.3e", x[i]);
        }
        printf("\n");
    }
    if (matrix_input(argc, n, k, r, s, tm, z, qb, a, b, argv, newcom) == -3) {
        free(a);
        free(ar);
        free(qb);
        if (r == 0) {
            free(b);
        }
        free(x);
        free(d);
        return -3;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double ans = 0;
    if (r == 0) {
        for (int i = 1; i < n; i += 2) {
            ans += x[i] * x[i];
        }
        ans = sqrt(ans);
    }
    residual(a, qb, b, x, n, r, s, tm, newcom);
    if (r == 0) {
        printf("Error: %10.3e\n", ans);
    }
    if (tm == 0) {
        MPI_Comm_free(&newcom);
    }
    MPI_Finalize();
    free(a);
    free(ar);
    if (r == 0) {
        free(b);
    }
    free(qb);
    free(d);
    free(x);
    return 0;
}
