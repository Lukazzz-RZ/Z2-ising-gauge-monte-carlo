#include "../include/ising_gauge.h"

// FUNCIONES RELACIONADAS CON LA ACTUALIZACIÓN DEL SISTEMA (MONTE CARLO)

// Cambia el signo del link en la dirección dada
inline void flip_link(int x, int y, int z, int dir) { // inline para optimizar
    if (dir == 0) Matriz[x][y][z].x *= -1;
    if (dir == 1) Matriz[x][y][z].y *= -1;
    if (dir == 2) Matriz[x][y][z].z *= -1;
}

// Calcula dE para un enlace dado
int calc_deltaE(int x, int y, int z, int d) {
    int DeltaE = 0;

    if (d == 0) {
        DeltaE = 2 * (Loops(1, x, y, z, 0) + Loops(1, x, y, z, 1) +
            Loops(1, x, y, z + off_n[z], 0) + Loops(1, x, y + off_n[y], z, 1));
    }
    if (d == 1) {
        DeltaE = 2 * (Loops(1, x, y, z, 1) + Loops(1, x, y, z, 2) +
            Loops(1, x + +off_n[x], y, z, 1) + Loops(1, x, y, z + off_n[z], 2));
    }
    if (d == 2) {
        DeltaE = 2 * (Loops(1, x, y, z, 2) + Loops(1, x, y, z, 0) +
            Loops(1, x + +off_n[x], y, z, 0) + Loops(1, x, y + off_n[y], z , 2));
    }

    return DeltaE;
}

// Un paso de Metropolis
int metropolis() {
    int n_x = (int)(L * N_Rand());
    int n_y = (int)(L * N_Rand());
    int n_z = (int)(L * N_Rand());
    int n_dir = (int)(3 * N_Rand());
    //printf("n_x: %d ,n_y: %d ,n_z: %d , dir: %d\t", n_x, n_y, n_z, n_dir);
    // calcular dE
    int Delta_E = calc_deltaE(n_x, n_y, n_z, n_dir);

    // Normalización al rango [0..4]
    int ind = Delta_E / 4 + 2;

    double r = N_Rand();
	//printf("r=%.5f, prob=%.5f, dE = %d\n", r, probs[ind], Delta_E);
    // criterio de Metropolis
	if (r < probs[ind]) {
        flip_link(n_x, n_y, n_z, n_dir);
		//printf("Aceptado dE = %d\n", Delta_E);
        return Delta_E;
    }

	return 0; // No se ha aceptado el cambio
}

int metropolis_con_seleccion(int n_x, int n_y, int n_z, int n_dir) {
    //printf("n_x: %d ,n_y: %d ,n_z: %d , dir: %d\t", n_x, n_y, n_z, n_dir);
    // calcular dE
    int Delta_E = calc_deltaE(n_x, n_y, n_z, n_dir);

    // Normalización al rango [0..4]
    int ind = Delta_E / 4 + 2;

    double r = N_Rand();
    //printf("r=%.5f, prob=%.5f, dE = %d\n", r, probs[ind], Delta_E);
    // criterio de Metropolis
    if (r < probs[ind]) {
        flip_link(n_x, n_y, n_z, n_dir);
        //printf("Aceptado dE = %d\n", Delta_E);
        return Delta_E;
    }

    return 0; // No se ha aceptado el cambio
}

// Termaliza el sistema bajo un criterio de variación dado
void termalize(double cond, int term_steps, int saveMat) {

    _mkdir("Results");
    _mkdir("Results/Comprobaciones");
    _mkdir("Results/MatC");

    // ABRIR FICHERO
    char filename_term[256];
    sprintf(filename_term, "Results/Comprobaciones/TermalizacionRed%d_beta%.2f.txt", L, beta);
    FILE* f = fopen(filename_term, "w");
    if (f == NULL) {
        printf("Error al abrir el fichero de termalización.\n");
        return;
    }

	// INICIO DE LA TERMALIZACIÓN
	int E0 = H();   // Energía inicial
    double avg_dE = cond + 1;
    int attempt = 0;

	// GUARDAR VALORES INICIALES
    fprintf(f, "%d\t%d\t%f\n", attempt * term_steps, E0, NAN);

	// ARRAY CIRCULAR PARA GUARDAR LOS ÚLTIMOS N_BLOCKS PROMEDIOS DE dE
    double last_dE[N_BLOCKS] = { 0 };
    int idx = 0;

	// BUCLE DE TERMALIZACIÓN (HASTA CONVERGENCIA)
    while (fabs(avg_dE) > cond && attempt < MAX_TERM_ATTP) {
        int dE_total = 0;

		// REALIZAR term_steps PASOS DE METROPOLIS
        for (int step = 0; step < term_steps; step++) {
            int dE = metropolis();
            dE_total += dE;
            E0 += dE;
        }

		// CALCULAR PROMEDIO DE dE EN LOS ÚLTIMOS term_steps PASOS (UN BLOQUE)
        double block_avg = (double)dE_total / term_steps;
        last_dE[idx % N_BLOCKS] = block_avg;
        idx++;

		// CALCULAR PROMEDIO DE LOS ÚLTIMOS N_BLOCKS BLOQUES
        int blocks_to_average = (idx < N_BLOCKS) ? idx : N_BLOCKS;
        double sum = 0;
        for (int i = 0; i < blocks_to_average; i++) sum += (last_dE[i]);
        avg_dE = sum / blocks_to_average;

		// SE ESCRIBE EN EL FICHERO EL PROMEDIO DE ENERGÍA DEL BLOQUE
        fprintf(f, "%d\t%d\t%.8f\n", attempt * term_steps, E0, block_avg);

		// GUARDAR CONFIGURACIÓN PARA VISUALIZACIÓN
        if (!saveMat) Mat_bin(attempt * term_steps);
        attempt++;
    }

	// CALCULAR extar NUEVOS BLOQUES AUNQUE YA HAYA CONVERGIDO (POR CLARIDAD)
    for (int extar = 0; extar < 50; extar++) {
        int dE_total = 0;
        for (int step = 0; step < term_steps; step++) {
            int dE = metropolis();
            dE_total += dE;
            E0 += dE;
        }
        double block_avg = (double)dE_total / term_steps;
		// GUARDAR VALORES EN EL FICHERO
        fprintf(f, "%d\t%d\t%.8f\n", (attempt + extar) * term_steps, E0, block_avg);
    }
    fclose(f);

	// CASO DE NO CONVERGENCIA/CONVERGENCIA
    if (attempt == MAX_TERM_ATTP)
        printf("Aviso: la termalizacion no ha convergido\n");
    else
        printf("Termalizacion completada en %d pasos, dE/N promedio = %.10f\n",
            attempt * term_steps, avg_dE);
}

// Cambia una configuracion por otra con misma energia
void overrelax() {
    int dE;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            for (int z = 0; z < L; z++){
                for (int d = 0; d < 3; d++) {
                    dE = calc_deltaE(x, y, z, d);
                    if (dE == 0) flip_link(x, y, z, d);
                }
            }
        }
    }
}
