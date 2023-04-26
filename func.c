#include "func.h"

double func(int k, int n, int i, int j)
{
    int max;
    double p;
    max = i;
    if (max < j) {
        max = j;
    }
    if (k == 1) {
        p = n - max + 1;
    } else if (k == 2) {
        p = max;
    } else if (k == 3) {
        p = abs(i - j);
    } else {
        p = 1.0 / (i + j - 1);
    }
    return p;
}
int matrix_input(int argc, int n, int k, int r, int s, int tm, int z, double *qb, double *a,
            double *b, char *argv[], MPI_Comm newcom)
{
    int p;
    if (argc == MAX) {
        for (int i = 0; i < n; i++) {
            if (r == 0) {
                for (int j = 0; j < n; j++) {
                    b[j] = func(k, n, i + 1, j + 1);
                }
                qb[i] = 0;
                for (int j = 0; j < (n + 1) / 2; j++) {
                    qb[i] += b[2 * j];
                }
            }
            if (s <= n) {
                for (int j = 0; j < n - z; j += s) {
                    MPI_Scatter(b + j, 1, MPI_DOUBLE,
                                a + i * (n / s + 1) + j / s, 1, MPI_DOUBLE, 0,
                                MPI_COMM_WORLD);
                }
            }
            if (tm == 0) {
                MPI_Scatter(b + n - z, 1, MPI_DOUBLE,
                            a + i * (n / s + 1) + n / s, 1, MPI_DOUBLE, 0,
                            newcom);
            }
        }
    } else {
        FILE *fp = NULL;
        int pp = 0;
        fp = fopen(argv[MAX], "r");
        if (!fp) {
            printf("указанного файла не существует\n");
            return -3;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (r == 0) {
                    p = fscanf(fp, "%lf", &b[j]);
                    if (p != 1) {
                        printf("не удалось считать матрицу из файла\n");
                        pp = 1;
                        break;
                    }
                }
            }
            if (r == 0) {
                qb[i] = 0;
                for (int j = 0; j < (n + 1) / 2; j++) {
                    qb[i] += b[2 * j];
                }
            }
            MPI_Bcast(&pp, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (pp) {
                fclose(fp);
                return -3;
            }
            if (s <= n) {
                for (int j = 0; j < n - z; j += s) {
                    MPI_Scatter(b + j, 1, MPI_DOUBLE,
                                a + i * (n / s + 1) + j / s, 1, MPI_DOUBLE, 0,
                                MPI_COMM_WORLD);
                }
            }
            if (tm == 0) {
                MPI_Scatter(b + n - z, 1, MPI_DOUBLE,
                            a + i * (n / s + 1) + n / s, 1, MPI_DOUBLE, 0,
                            newcom);
            }
        }
        fclose(fp);
    }
    return 0;
}
struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec - start.tv_nsec) < 0) {
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = GIGA_MODIFIER + end.tv_nsec - start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return temp;
}

int Cholesky(double *a, double *ar, double *q, double *x, double *d, double *b, int n, int r, int s,
       int tm, MPI_Comm newcom)
{
    struct timespec time_proc_start;
    struct timespec time_proc_end;
    clock_gettime(CLOCK_MONOTONIC, &time_proc_start);
    int u = n / s + 1;
    MPI_Status status;
    int k;
    int j;
    int i;
    int t;
    int z1;
    int m;
    double p;
    double tmp;
    for (i = 0; i < n; i++) {
        tmp = 0;
        t = 0;
        m = i;
        if (r == i % s) {
            tmp = fabs(a[i * u + i / s]);
            for (j = i + 1; j < n; j++) {
                if (fabs(a[j * u  + i / s]) > tmp) {
                    tmp = fabs(a[j * u  + i / s]);
                    m = j;
                }
            }
        }
        MPI_Bcast(&m, 1, MPI_INT, i % s, MPI_COMM_WORLD);
        z1 = 0;
        if (r < i % s) {
            z1 = 1;
        }
        if (m != i) {
            for (j = i / s + z1; j < u - tm; j++) {
                p = a[i * u + j];
                a[i * u + j] = a[m * u + j];
                a[m * u + j] = p;
            }
            p = q[i];
            q[i] = q[m];
            q[m] = p;
        }
        if (r == i % s) {
            if (fabs(a[i * u + i / s]) < 1e-50) {
                t = 1;
            }
        }
        MPI_Bcast(&t, 1, MPI_INT, i % s, MPI_COMM_WORLD);
        if (t == 1) {
            return -4;
        }
        for (j = 0; j < n; j++) {
            if (j != i) {
                if (r == i % s) {
                    p = a[j * u + i / s] / a[i * u + i / s];
                }
                MPI_Bcast(&p, 1, MPI_DOUBLE, i % s, MPI_COMM_WORLD);
                for (k = i / s + z1; k < u - tm; k++) {
                    a[j * u + k] -= a[i * u + k] * p;
                }
                q[j] -= q[i] * p;
            }
        }
    }

    for (i = 0; i < n; i++) {
        tmp = 0;
        if (r == i % s) {
            ar[i * u + i / s] = a[i * u + i / s];
            for (j = 0; j < i; j++) {
                ar[i * u + i / s] -= a[j * u + i / s] * a[j * u + i / s] * d[j];
            }
            tmp = fabs(ar[i * u + i / s]);
            p = ar[i * u + i / s] / tmp;
            ar[i * u + i / s] = sqrt(tmp);
        }
        MPI_Bcast(&p, 1, MPI_DOUBLE, i % s, MPI_COMM_WORLD);
        d[i] = p;
        z1 = 0;
        if (r < (i + 1) % s) {
            z1 = 1;
        }
        if (r == i % s) {
            for (k = 0; k < i + 1; k++) {
                x[k] = ar[k * u + i / s];
            }
        }
        MPI_Bcast(x, n, MPI_DOUBLE, i % s, MPI_COMM_WORLD);
        for (j = (i + 1) / s + z1; j < u - tm; j++) {
            tmp = a[i * u + j];
            for (k = 0; k < i; k++) {
                tmp -= x[k] * d[k] * ar[k * u + j];
            }
            ar[i * u + j] = tmp / (x[i] * d[i]);
        }
    }
    for (i = 0; i < n; i++) {
        if (r == i % s) {
            tmp = q[i];
            for (j = 0; j < i; j++) {
                tmp -= ar[j * u + i / s] * x[j];
            }
            p = tmp / ar[i * u + i / s];
        }
        MPI_Bcast(&p, 1, MPI_DOUBLE, i % s, MPI_COMM_WORLD);
        x[i] = p;
    }
    int z = n % s;
    for (k = 1; k <= n; k++) {
        if (s <= n) {
            for (j = 0; j < n - z; j += s) {
                MPI_Allgather(ar + (n - k) * u + j / s, 1, MPI_DOUBLE, q + j, 1,
                              MPI_DOUBLE, MPI_COMM_WORLD);
            }
            if (z != 0) {
                if (tm == 0) {
                    MPI_Gather(ar + (n - k) * u + n / s, 1, MPI_DOUBLE, b, 1,
                               MPI_DOUBLE, 0, newcom);
                }
                if (r == 0) {
                    for (j = 1; j < s; j++) {
                        MPI_Send(b, z, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
                    }
                    for (j = z; j >= 1; j--) {
                        q[n - j] = b[z - j];
                    }
                } else {
                    MPI_Recv(q + n - z, z, MPI_DOUBLE, 0, MPI_ANY_TAG,
                             MPI_COMM_WORLD, &status);
                }
            }
        } else {
            if (tm == 0) {
                MPI_Allgather(ar + (n - k) * u, 1, MPI_DOUBLE, q, 1, MPI_DOUBLE,
                              newcom);
            }
        }
        if (r == 0) {
            tmp = x[n - k];
            for (j = n - k + 1; j < n; j++) {
                tmp -= q[j] * x[j];
            }
            x[n - k] = tmp / (q[n - k] * d[n - k]);
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &time_proc_end);
    time_proc_end = diff(time_proc_start, time_proc_end);
    printf(
        "Time of process %d: %lf s\n", r,
        (double)(time_proc_end.tv_sec + time_proc_end.tv_nsec * NANO_MODIFIER));
    return 0;
}
void residual(const double *a, const double *q, double *b, double *x, int n, int r,
         int s, int tm, MPI_Comm newcom)
{
    int u = n / s + 1;
    int z = n % s;
    double norm = 0.0;
    double tmp = 0.0;
    for (int j = 0; j < n; j++) {
        if (s <= n) {
            for (int i = 0; i < n - z; i += s) {
                MPI_Gather(a + j * u + i / s, 1, MPI_DOUBLE, b + i, 1,
                           MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            if (tm == 0) {
                MPI_Gather(a + j * u + n / s, 1, MPI_DOUBLE, b + n - z, 1,
                           MPI_DOUBLE, 0, newcom);
            }
        } else {
            if (tm == 0) {
                MPI_Gather(a + j * u, 1, MPI_DOUBLE, b, 1, MPI_DOUBLE, 0,
                           newcom);
            }
        }
        if (r == 0) {
            norm += q[j] * q[j];
            double t = 0;
            for (int i = 0; i < n; i++) {
                t += b[i] * x[i];
            }
            tmp += (q[j] - t) * (q[j] - t);
            t = 0.0;
        }
    }
    if (r == 0) {
        printf("Residual: \n%10.3e\n", sqrt(tmp / norm));
    }
}

void out(double *a, double *b, double *qb, int n, int m, int r, int s, int tm,
         MPI_Comm newcom)
{
    int z = n % s;
    for (int j = 0; j < n; j++) {
        if (s <= n) {
            for (int i = 0; i < n - z; i += s) {
                MPI_Gather(a + j * (n / s + 1) + i / s, 1, MPI_DOUBLE, b + i, 1,
                           MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            if (tm == 0) {
                MPI_Gather(a + j * (n / s + 1) + n / s, 1, MPI_DOUBLE,
                           b + n - z, 1, MPI_DOUBLE, 0, newcom);
            }
        } else {
            if (tm == 0) {
                MPI_Gather(a + j * (n / s + 1), 1, MPI_DOUBLE, b, 1, MPI_DOUBLE,
                           0, newcom);
            }
        }
        if (r == 0 && j < m) {
            for (int i = 0; i < m; i++) {
                printf(" %10.3e", b[i]);
            }
            printf(" %10.3e", qb[j]);
            printf("\n");
        }
    }
    if (r == 0) {
        printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

int input(char *t)
{
    char *end;
    int b = (int)strtol(t, &end, DECIMAL);
    if (*end != '\0') {
        printf("не удалось считать число\n");
        return -1;
    }
    return b;
}
