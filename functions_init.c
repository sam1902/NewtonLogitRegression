//
// Created by Samuel Prevost on 23/12/2019.
//

#include "functions_init.h"

void compute_hessian(double** hess, const double** A, const double** B, const double* theta){
    for (int i = 0; i < N_DIM; ++i) {
        for (int j = 0; j < N_DIM; ++j) {
            hess[i][j] = dFlk(A, B, theta, i, j);
        }
    }
}

void compute_grad(double* grad, const double** A, const double** B, const double* theta){
    for (int i = 0; i < N_DIM; ++i) {
        grad[i] = dFl(A, B, theta, i);
    }
}

double F(const double** A, const double** B, const double* T){
    double res = 0;
    for (int i = 0; i < DATA_SIZE; ++i) {
        res += log(1+exp(-prod(T, A[i]))) + log(1+exp(prod(T, B[i])));
    }
    return res;
}

double dFl(const double **A, const double **B, const double *T, int l){
    double res = 0;
    for (int i = 0; i < DATA_SIZE; ++i) {
        res += A[i][l]*(sig(prod(T, A[i])) - 1) + B[i][l]*sig(prod(T, B[i]));
    }
    return res;
}

double dFlk(const double** A, const double** B, const double* T, int l, int k){
    double res = 0;
    for (int i = 0; i < DATA_SIZE; ++i) {
        double sig_prod_theta_a_i = sig(prod(T, A[i]));
        double sig_prod_theta_b_i = sig(prod(T, B[i]));
        res += A[i][l] * A[i][k] * sig_prod_theta_a_i * (1 - sig_prod_theta_a_i);
        res += B[i][l] * B[i][k] * sig_prod_theta_b_i * (1 - sig_prod_theta_b_i);
    }
    return res;
}