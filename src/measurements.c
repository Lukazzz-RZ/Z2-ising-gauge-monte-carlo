#include "../include/ising_gauge.h"

// Calcula productos de links en un loop nxn
int Loops(int n, int x, int y, int z, int d) {
    if (n >= L) return 0; 

    int loop = 1;
 
    int v[3] = {x ,y ,z};
    int d1, d2;
    if (d == 1) { d1 = 0; d2 = 1; }      // plano XY
    else if (d == 2) { d1 = 1; d2 = 2; } // plano YZ
    else { d1 = 0; d2 = 2; }             // plano ZX

    // borde 1: +d1
    for (int i = 0; i < n; i++) {
        loop *= get_link(v[0], v[1], v[2], d1);
        v[d1] = (v[d1] + off_p[v[d1]] + L) % L;
    }

    // borde 2: +d2
    for (int i = 0; i < n; i++) {
        loop *= get_link(v[0], v[1], v[2], d2);
        v[d2] = (v[d2] + off_p[v[d2]] + L) % L;
    }

    // borde 3: -d1
    for (int i = 0; i < n; i++) {
        v[d1] = (v[d1] + off_n[v[d1]] + L) % L;
        loop *= get_link(v[0], v[1], v[2], d1);
    }

    // borde 4: -d2
    for (int i = 0; i < n; i++) {
        v[d2] = (v[d2] + off_n[v[d2]] + L) % L;
        loop *= get_link(v[0], v[1], v[2], d2);
    }

    return loop;
}

// Calcula loops desde nMin hasta nMax
void LoopsHastaNmax(int nMin, int nMax, int* ELoops, int x, int y, int z, int d) {
    for (int i = 0; i <= nMax - nMin; i++) {
        int n = nMin + i;
        ELoops[i] = Loops(n, x, y, z, d);
    }
}

// Calcula la energa de la red
long int H(void) {
    int H0 = 0;

    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            for (int z = 0; z < L; z++) {
                for (int i = 0; i < 3; i++) {
                    H0 += Loops(1,x, y, z, i);
                }
            }
        }
    }
    return -H0;
}