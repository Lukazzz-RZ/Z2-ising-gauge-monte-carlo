#include "../include/ising_gauge.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main() {

    ////// INICIALIZACIÓN //////
    Ini_N_Rand(time(NULL));
    precalc_exp(probs);
    printf("PROBABILIDADES DE ACEPTACION A á = %f:\n\n", beta);
    for (int i = 0; i < 5; i++) {
        printf("prob de dE = %d: %f\n", (-4 + 4 * (i - 1)), min(1,probs[i]));
    }
    precalc_offset(off_p, off_n);
    ini_Matriz_rand();
    //ini_Matriz_const(1);

    int loop_min = 1;
    int loop_max = 10;
    int n_loop_sizes = loop_max - loop_min + 1;
    long int* loop_values = malloc(sizeof(int) * n_loop_sizes);
    long double* global_sum = calloc(n_loop_sizes, sizeof(double));
    long int* loop_values_square = malloc(sizeof(int) * n_loop_sizes);
    long double* global_sum_square = calloc(n_loop_sizes, sizeof(double));

    int n_configs = 1000;
    int total_links = L * L * L * 3;
    int steps_per_config = 200000;

    double termo_eps = 1e-24;
    int termo_bloque = 2500;                 
    

    ////// ARCHIVOS DE SALIDA //////
    _mkdir("Results");
    _mkdir("Results/Comprobaciones");

    char filename[256];
    char filename2[256];
    char filename3[256];
    sprintf(filename, "Results/Comprobaciones/LoopsRed%d_beta%.2f.txt", L, beta);
    sprintf(filename2, "Results/Comprobaciones/Energia%d_beta%.2f.txt", L, beta);
    sprintf(filename3, "Results/Comprobaciones/TermalizacionRed%d_beta%.2f.txt", L, beta);

    FILE* fout = fopen(filename, "w");
    FILE* fout2 = fopen(filename2, "w");
    if (!fout || !fout2) {
        fprintf(stderr, "Error al abrir archivos de salida.\n");
        return EXIT_FAILURE;
    }

    fprintf(fout, "# Loops de Wilson promedio\n# L = %d, beta = %.3f, n_configs = %d\n# loop_size\t<W(loop)>\n",
        L, beta, n_configs);
    fprintf(fout2, "# Energia de configuraciones medidas\n# L = %d, beta = %.3f, n_configs = %d\n# E\t<E>\n",
        L, beta, n_configs);

    ////// TERMALIZACIÓN //////
    printf("\nCOMENZANDO TERMALIZADO: \n\n");
    termalize(termo_eps, termo_bloque, 1);
    overrelax();
    int E0 = H();
    printf("\nCOMENZANDO ADQUISICION: \n\n");
    ////// BUCLE PRINCIPAL //////
    int interval = (n_configs >= 20) ? n_configs / 20 : 1;
    int dE = 0;
    long double avg_E = E0;
    fprintf(fout2, "%f\t%.6f\n", (double)E0, (double)avg_E);

    int N = 0;

    for (int config = 0; config < n_configs; config++) {
        if (config % 5 == 0) overrelax();
        for (int step = 0; step < steps_per_config; step++) {
            dE = metropolis();
            // Seguimiento Energias
            E0 += dE;
            avg_E = ((config * steps_per_config + step) * avg_E + E0) /
                (config * steps_per_config + step + 1);
            fprintf(fout2, "%f\t%.6f\n", (double)E0, (double)avg_E);
        }

        // ALGUNOS LOOPS DE UNA CONFIGURACION (AL AZAR XD)
		int loops_per_config = L*L;
        for (int i = 0; i < loops_per_config; i++) {
            int n_x = (int)(L * N_Rand());
            int n_y = (int)(L * N_Rand());
            int n_z = (int)(L * N_Rand());
            int n_dir = (int)(3 * N_Rand());
			LoopsHastaNmax(loop_min, loop_max, loop_values, n_x, n_y, n_z, n_dir);
            for (int i = 0; i < n_loop_sizes; i++) {
                loop_values_square[i] = loop_values[i] * loop_values[i];
            }
            for (int i = 0; i < n_loop_sizes; i++) {
                global_sum[i] += loop_values[i];
                global_sum_square[i] += loop_values_square[i];
            }
			N++;
        }

		// TODOS LOS LOOPS DE UNA CONFIGURACION
     /*   for (int x = 0; x < L; x += 1) {
            for (int y = 0; y < L; y += 1) {
                for (int z = 0; z < L; z += 1) {
                    for (int dir = 0; dir < 3; dir++) {
                        LoopsHastaNmax(loop_min, loop_max, loop_values, x, y, z, dir);
                        for (int i = 0; i < n_loop_sizes; i++) {
                            global_sum[i] += loop_values[i];
                            
                        }
                        N++;
                    }
                }
            }
        }*/

        if ((config % interval) == 0) {
            printf("Progreso: %d%%\n", (config * 100) / n_configs);
        }
    }


    //////
    for (int i = 0; i < n_loop_sizes; i++) {
        double wloop_avg = global_sum[i] / (double)(N);
        double wloop_square_avg = global_sum_square[i] / (double)(N);
        fprintf(fout, "%d\t%.6f\t%.6f\n", loop_min + i, wloop_avg,sqrt(wloop_square_avg- wloop_avg* wloop_avg)/sqrt((double)N));
    }

    fclose(fout);
    fclose(fout2);
    free(loop_values);
    free(global_sum);

    printf("\nSimulacion completada.\n");
    printf("Resultados guardados en:\n  %s\n  %s\n", filename, filename2);

    ///// EJECUCIÓN DE SCRIPTS PYTHON /////
    char comando1[256];
    char comando2[256];

    // Pasar archivo de energía a Check_termalizado.py
    snprintf(comando1, sizeof(comando1),
        "python3 Results/Check_termalizado.py \"%s\"", filename3);

    if (beta <= 0.7613) {
        snprintf(comando2, sizeof(comando2),
            "python3 Results/Check_arealaw.py \"%s\"", filename);
    }
    else {
        snprintf(comando2, sizeof(comando2),
            "python3 Results/Check_perilaw.py \"%s\"", filename);
    }
    printf("\nEjecutando analisis:\n  %s\n  %s\n", comando1, comando2);

    int status1 = system(comando1);
    int status2 = system(comando2);

    if (status1 != 0)
        fprintf(stderr, "[!] El script Check_termalizado.py devolvio un error.\n");
    if (status2 != 0)
        fprintf(stderr, "[!] El script Check_arealaw.py devolvio un error.\n");

    printf("\nAnalisis completado.\n");
    return 0;
}
