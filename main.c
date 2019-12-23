#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

/* Toutes les directives de préprocesseur définissant
 * des constantes sont déclarées dans macros.h */
#include "marcros.h"

#include "gauss.h" /* Contient la méthode gauss */
#include "vector_utils.h" /* Contient les fonctions de manipulation de vecteurs */
#include "functions_init.h" /* Contient la fonction objectif, ses dérivées et des aides à leurs calculs */

double** load_data_file(const char* file_name);

int main(){
    const double** A = (const double**) load_data_file("Base_apprentissage_Notes_Admis.dat");
    const double** B = (const double**) load_data_file("Base_apprentissage_Notes_NonAdmis.dat");

    double *C0 = init_vect(N_DIM);
    double* ones_ = init_vect(N_DIM);
    double* theta = init_vect(N_DIM);
    for (int i = 0; i < N_DIM; ++i) {
        theta[i] = THETA_START;
        ones_[i] = 1;
    }
    const double* ONES = ones_;

    // Dans le cas où la hessienne n'est plus inversible, réessaye en repartant de -theta
    bool tried_other_side = false;

    double* grad = init_vect(N_DIM);
    double** hess = init_mat(N_DIM, N_DIM);

    int iteration = 0;
    do {
        for (int i = 0; i < N_DIM; ++i) {
            C0[i] = 0;
            grad[i] = 0;
            for (int j = 0; j < N_DIM; ++j) {
                hess[i][j] = 0;
            }
        }
        compute_hessian(hess, A, B, theta);
        compute_grad(grad, A, B, theta);

        for (int i = 0; i < N_DIM; ++i) {
            C0[i] = grad[i] - prod(hess[i], theta);
        }

        if(!gauss(N_DIM, hess, C0, theta)){
            printf("#################################################\n");
            if (!tried_other_side) {
                printf("# Echec ! Nouvel essai en partant de -theta/|theta|*0.01...\n");
                tried_other_side = true;
                int _discard = 0;
                for (int i = 0; i < N_DIM; ++i) {
                    theta[i] = -frexp(theta[i], &_discard)*0.01;
                }
                printf("# Nouveau theta: \n# ");
                print_vect_precise(theta);
            } else {
                printf("# Echec total..\n");
                exit(1);
            }
            printf("#################################################\n");
        }
        iteration++;
        if (iteration & 1){  // Print que les itérations paires
            printf("########### it: %d ###########\n", iteration);
            printf("\tHess:\n");
            print_mat(hess);
            printf("\tGrad: \t");
            print_vect(grad);
            printf("\tGrad sum: %.5f\n", fabs(prod(grad, ONES)));
            printf("\tC0:\t");
            print_vect(C0);
            printf("\tTheta:\t");
            print_vect(theta);
        }
        // Stop quand les gradients sont nuls
    } while (fabs(prod(grad, ONES)) > EPS );

    printf("Convergence à %.15f près en %d iterations\n", EPS, iteration);
    printf("Theta optimal:\n");
    print_vect_precise(theta);
    return 0;
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cert-err34-c"
/***
 * Charges le fichier de données `file_name` en mémoire dans une matrice avec une colonne de 1 en plus.
 * @param file_name : le fichier de données
 * @return une matrice avec une colonne de 1 en plus
 */
double** load_data_file(const char* file_name){
    double** A;
    A=(double**) malloc(sizeof(double*) * DATA_SIZE);
    for (int i=0; i < DATA_SIZE; i++){
        A[i]=(double*) malloc(sizeof(double) * N_DIM);
    }

    FILE *f1=fopen(file_name,"r");
    if (!f1){printf("Problème avec le fichier %s\n", file_name); return 0;}
    for (int i=0; i < DATA_SIZE; i++) {
        for (int k = 0; k < N_DIM - 1; k++) {
            fscanf(f1, "%lf", &(A[i][k]));
        }
        A[i][N_DIM - 1] = 1;  /* Permet de creer une constante dans le produit scalaire theta^T.A */
    }
    fclose(f1);
    return A;
}
#pragma clang diagnostic pop