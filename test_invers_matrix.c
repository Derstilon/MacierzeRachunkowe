#include <stdio.h>
// #include <malloc.h>
#include <math.h>
#include <stdlib.h>

#include "parameters.h"

void MY_MMult(double **, double **, double **, int, int, long long *, long long *);
void copy_matrix(int, int, double **, int, double **, int);
void random_matrix(int, int, double **, int);
double compare_matrices(int, int, double **, int, double **, int);
void print_matrix(int m, int n, double **a, int lda);
double dclock();
void inverse_matrix(double **, double **, int, long long *, long long *);

double fabs(double x) {
    if (x > 0) {
        return x;
    } else {
        return -x;
    }
}

int main() {
    int l, i, j;
    double **a, **a_inv, **b;
    double start_time, end_time, best_time;
    long long mult_count, add_count, tmp_mult_count, tmp_add_count;
    int status = 1;

    for (l = 2; l <= pow(2, LLAST); l *= 2) {
        a = (double**) malloc(l * sizeof(double*));
        a_inv = (double**) malloc(l * sizeof(double *));
        b = (double **)malloc(l * sizeof(double *));
        for(i = 0;i<l;i++){
            a[i] = (double*) malloc(l * sizeof(double));
            a_inv[i] = (double *) malloc(l * sizeof(double));
            b[i] = (double *)malloc(l * sizeof(double));
        }

        random_matrix(l, l, a, 0);

        for (i = 0; i < NREPEATS; i++){
            start_time = dclock();
            inverse_matrix(a, a_inv, l, &mult_count, &add_count);
            end_time = dclock();
            if (i == 0) {
                best_time = end_time - start_time;
            } else {
                best_time = (end_time - start_time < best_time ? end_time - start_time : best_time);
            }
        }

        /// Error check
        //print_matrix(l, l, a, l);
        //print_matrix(l, l, a_inv, l);

        MY_MMult(a, a_inv, b, l, 8, &tmp_mult_count, &tmp_add_count);
        //print_matrix(l, l, b, l);
        status = 1;
        for (i = 0; i < l; i++) {
            for (j = 0; j < l; j++) {
                if (i == j) {
                    if (fabs(1 - b[i][j]) > EPS) status = 0;
                } else {
                    if (fabs(b[i][j]) > EPS) status = 0;
                }
            }
        }

        printf("%d %d %le %lld %lld\n", l, status, best_time, mult_count, add_count);
        fflush(stdout);

        for (i = 0; i < l; i++) {
            free(a[i]);
            free(a_inv[i]);
            free(b[i]);
        }
        free(a);
        free(a_inv);
        free(b);
    }
    return 0;
}