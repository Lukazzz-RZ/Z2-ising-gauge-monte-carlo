#include "../include/ising_gauge.h"

// FUNCIONES RELACIONADAS CON LA INICIALIZACIN DEL SISTEMA Y PRECALCULOS


//Inicializado de variables globales
Nodo Matriz[L][L][L];
double beta = 0.8;        // Inversa de la temperatura
int MAX_TERM_ATTP = 100000;  // Nmero mximo de intentos de termalizacin
long double probs[5];      // Array para almacenar los valores precalculados de exp(-beta*dE)
int off_p[L];              // Array para almacenar los offsets peridicos positivos
int off_n[L];              // Array para almacenar los offsets peridicos negativos


// Precalcula valores de exp(-beta*dE) para dE = -8, -4, 0, 4, 8
void precalc_exp(double* probs) {
    for (int i = 0; i < 5; i++) {
        probs[i] = exp(-beta * (-8 + 4 * i));
    }
    return;
}

// Precalcula offsets peridicos para una red de tamao L
void precalc_offset(int* off_p, int* off_n) {
    for (int i = 0; i < L - 1; i++) {
        off_p[i] = 1;
        off_n[i + 1] = -1;
    }
    off_p[L - 1] = -(L - 1);
    off_n[0] = (L - 1);
    return;
}

// Asigna valores aleatorios (+1 o -1) a todos los links de la matriz
void ini_Matriz_rand() {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                Matriz[i][j][k].x = 2 * ((int)(N_Rand() + 0.5)) - 1;
                Matriz[i][j][k].y = 2 * ((int)(N_Rand() + 0.5)) - 1;
                Matriz[i][j][k].z = 2 * ((int)(N_Rand() + 0.5)) - 1;
            }
        }
    }
}

void ini_Matriz_const(int val) {
    if (val != 1 && val != -1) {
        printf("Error: ini_Matriz_const solo acepta 1 o -1\n");
        return;
    }

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                Matriz[i][j][k].x = val;
                Matriz[i][j][k].y = val;
                Matriz[i][j][k].z = val;
            }
        }
    }
}