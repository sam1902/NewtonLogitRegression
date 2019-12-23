//
// Created by Samuel Prevost on 23/12/2019.
//

#ifndef REGLOG_VECTOR_UTILS_H
#define REGLOG_VECTOR_UTILS_H

#include "marcros.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void print_mat(double** mat);
void print_vect(double* vect);
void print_vect_precise(double* vect);
double* init_vect(int n);
double** init_mat(int n, int m);
double sig(double z);
double prod(const double* T, const double* X);
double norm2(const double* X1, const double* X2);

#endif //REGLOG_VECTOR_UTILS_H
