#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <direct.h> 
#include <time.h>

#define NormRANu (2.3283063671E-10F)
#define L 24
#define N_BLOCKS 100

typedef struct {
    int x;
    int y;
    int z;
} Nodo;

/// Paris y Rapuano
unsigned int irr[256];
unsigned int ir1;
unsigned char ind_ran, ig1, ig2, ig3;


/// Variables globales del modelo
extern double beta;                     // Inversa de la temperatura
extern int MC_STEPS;                    // Número de pasos de Monte Carlo
extern int MAX_TERM_ATTP; 			    // Número máximo de intentos de termalización
extern Nodo Matriz[L][L][L];            // Matriz que contiene los nodos del sistema
extern long double probs[5];            // Array para almacenar los valores precalculados de exp(-beta*dE)
extern int off_p[L];                    // Array para almacenar los offsets periódicos positivos
extern int off_n[L];                    // Array para almacenar los offsets periódicos negativos



double N_Rand(void);                    
void Ini_N_Rand(int SEMILLA);
void Mat_bin(int n);
int get_link(int x, int y, int z, int dir);

void precalc_exp(double* probs);
void precalc_offset(int* off_p, int* off_n);
void ini_Matriz_rand();
void ini_Matriz_const(int val);

int metropolis();
void termalize(double cond, int term_steps, int saveMat);
void overrelax();

long int H(void);
int Loops(int n, int x, int y, int z, int d);