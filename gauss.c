//
// Created by Samuel Prevost on 23/12/2019.
//

#include "gauss.h"


/***
 * Resolution du système
 * C.x + C0 = 0
 * @param n_dim dimension count
 * @param C matrix of shape (N_DIM, N_DIM)
 * @param C0 vector of shape (N_DIM)
 * @param x output vector, contains solution if it exists
 * @return 1 si le systeme a une unique solution et 0 sinon
 */
int gauss(int n_dim, double** C, double* C0, double* x){
    int imin;
    double pivot;
    double sum, valmin, temp1, temp2;

    for(int k=0; k < n_dim - 1; k++) {
        /* Recherche de l'element de valeur absolue minimum non nulle
           dans la colonne k et d'indice i>=k.*/

        valmin = C[k][k]; imin = k;
        for(int i = k+1; i < n_dim; i++){
            if (valmin != 0){
                if (fabs(C[i][k]) < fabs(valmin) && C[i][k] != 0){
                    valmin = C[i][k];
                    imin = i;
                }
            } else {
                valmin = C[i][k];
                imin = i;
            }
        }

        /* Si l'element minimum est nul, on peut en déduire
           que le systeme n'a pas de solution ou n'a pas de solution
           unique.*/

        if (valmin == 0.) return 0;

        /* Sinon inversion des elements de la ligne imax avec les elements
           de la ligne k. On fait de meme avec le vecteur C0. */

        for(int j=0; j < n_dim; j++){
            temp1 = C[imin][j];
            C[imin][j] = C[k][j];
            C[k][j] = temp1;
        }

        temp2 = C0[imin];
        C0[imin] = C0[k];
        C0[k] = temp2;

        /* Reduction de la matrice par la methode d'élimination de Gauss */

        for(int i=k+1; i < n_dim; i++) {
            pivot = C[i][k]/C[k][k];
            for(int j=0; j < n_dim; j++) {
                C[i][j] = C[i][j] - pivot * C[k][j];
            }
            C0[i] = C0[i] - pivot*C0[k];
        }
    }

    if (C[n_dim - 1][n_dim - 1] == 0){
        return 0;
    }

    /* Deduction de l'unique solution */

    x[n_dim - 1] = -C0[n_dim - 1] / C[n_dim - 1][n_dim - 1];

    for(int i= n_dim - 2; i > -1; i--){
        sum = 0;
        for(int j= n_dim - 1; j > i; j--) {
            sum = sum + C[i][j] * x[j];
        }
        x[i] = (-C0[i] - sum)/C[i][i];
    }
    return 1;
}