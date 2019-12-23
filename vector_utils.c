//
// Created by Samuel Prevost on 23/12/2019.
//

#include "vector_utils.h"

void print_mat(double** mat){
    for (int i = 0; i < N_DIM; ++i) {
        for (int j = 0; j < N_DIM; ++j) {
            printf("%.5f, ", mat[i][j]);
        }
        printf("\n");
    }
}

void print_vect(double* vect){
    printf("[ ");
    for (int i = 0; i < N_DIM; ++i) {
        printf("%.3f, ", vect[i]);
    }
    printf(" ]\n");
}

void print_vect_precise(double* vect){
    printf("[ ");
    for (int i = 0; i < N_DIM; ++i) {
        printf("%.15f, ", vect[i]);
    }
    printf(" ]\n");
}

double* init_vect(int n){
    double* vect = (double*) malloc(sizeof(double) * n);
    return vect;
}

double** init_mat(int n, int m){
    double** mat = (double**) malloc(sizeof(double*) * n);
    for (int i = 0; i < m; i++) {
        mat[i] = init_vect(m);
    }
    return mat;
}

double sig(double z){
    return 1/(1+exp(-z));
}

double prod(const double* T, const double* X){
    double res = 0;
    for (int i = 0; i < N_DIM; ++i) {
        res += T[i] * X[i];
    }
    return res;
}

double norm2(const double* X1, const double* X2){
    double res = 0;
    for (int i = 0; i < N_DIM; ++i) {
        res += pow(X1[i] - X2[i], 2);
    }
    return sqrt(res);
}