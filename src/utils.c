#include "../include/ising_gauge.h"

// FUNCIONES MATEMTICAS AUXILIARES
double N_Rand(void){

    double r;

    ig1 = ind_ran-24;
    ig2 = ind_ran-55;
    ig3 = ind_ran-61;

    irr[ind_ran] = irr[ig1] + irr[ig2];
    ir1 = (irr[ind_ran]^irr[ig3]);
    ind_ran++;
    r=ir1*NormRANu;

    return r;
}

void Ini_N_Rand(int Seed){

    int INI,FACTOR,SUM,i;

    INI = Seed;
    FACTOR = 67397;
    SUM = 7364893;
    srand(Seed);

    for(i=0; i<256; i++){
        INI = (INI * FACTOR + SUM);
        irr[i] = INI;
    }
    ind_ran = ig1 = ig2 = ig3 = 0;

    return;
}

void Mat_bin(int n) {

    char nombre[35];
    snprintf(nombre, sizeof(nombre), "Results\\MatC\\conf%010d.bin", n);
    FILE* f;
    f = fopen(nombre, "wb");
    if (f == NULL) {
        printf("Error al abrir el fichero conf%010d.bin\n", n);
        return;
    }

    fwrite(Matriz, sizeof(Nodo), L * L * L, f);
    fclose(f);
}

int get_link(int x, int y, int z, int dir) { // inline para optimizar
    if (dir == 0) return Matriz[x][y][z].x;
    if (dir == 1) return Matriz[x][y][z].y;
    return Matriz[x][y][z].z;
}