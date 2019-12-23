//
// Created by Samuel Prevost on 23/12/2019.
//

#ifndef REGLOG_FUNCTIONS_INIT_H
#define REGLOG_FUNCTIONS_INIT_H

#include "marcros.h"
#include "vector_utils.h"
#include <math.h>

void compute_hessian(double** hess, const double** A, const double** B, const double* theta);
void compute_grad(double* grad, const double** A, const double** B, const double* theta);
double F(const double** A, const double** B, const double* T);
double dFl(const double **A, const double **B, const double *T, int l);
double dFlk(const double** A, const double** B, const double* T, int l, int k);

#endif //REGLOG_FUNCTIONS_INIT_H
