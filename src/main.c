#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "rk.h"
#include "pendulo.h"

// Protótipos para as novas funções de análise
void run_comparative_analysis();
void run_10_period_analysis();
void find_max_angle_for_error_threshold();
void generate_plot_data();

int main() {
    printf("Executando analise comparativa...\n");
    run_comparative_analysis();

    printf("Executando Análise de Desempenho para 10 Períodos...\n");
    run_10_period_analysis();
    
    printf("\nProcurando ângulo para erro < 0.001...\n");
    find_max_angle_for_error_threshold();

    printf("\nGerando arquivos de dados para gráficos...\n");
    generate_plot_data();
    printf("Arquivos de dados gerados.\n");

    return 0;
}


/**
Fa¸ca um comparativo do valor do per´ıodo calculado e do n´umero de passos, para diferentes
ˆangulos iniciais θ0:
*/
void run_comparative_analysis(void) {
    FILE *fp = fopen("output/analise_completa.csv", "w");

    if (fp == NULL) {
        perror("Erro ao abrir o arquivo de saída");
        return;
    }

    // Lista de ângulos iniciais para testar (em radianos)
    double theta0_vals[] = {0.1, 0.5, 1.0, 2.0, 3.0};
    int n_thetas = sizeof(theta0_vals) / sizeof(theta0_vals[0]);

    // Parâmetros para os métodos numéricos
    double h_vals[] = {0.01, 0.001, 0.0001};
    int n_h = sizeof(h_vals) / sizeof(h_vals[0]);
    double tol_adapt = 1e-7;
    double h0_adapt = 0.01;

    // Escreve o cabeçalho no arquivo
    fprintf(fp, "theta0,method,h,period,steps,error_vs_adapt\n");

    double T_analytic = analytic_period();

    // Loop principal sobre cada ângulo inicial
    for (int i = 0; i < n_thetas; ++i) {
        double theta0 = theta0_vals[i];
        double T_adaptive;
        int steps_adaptive;

        // 1. Solução Analítica Simplificada
        fprintf(fp, "%.2f,analytic,N/A,%.8f,0,N/A\n", theta0, T_analytic);

        // 2. Solução com Passo Adaptativo (referência de precisão)
        detect_period_adaptive(theta0, tol_adapt, h0_adapt, 1, &T_adaptive, &steps_adaptive, NULL);
        fprintf(fp, "%.2f,adaptive,%.1e,%.8f,%d,0.0\n", theta0, tol_adapt, T_adaptive, steps_adaptive);
        
        // 3. Soluções com Passo Constante
        for (int j = 0; j < n_h; ++j) {
            double h = h_vals[j];
            int steps_const;
            double T_const = detect_period_constant(theta0, h, 1, &steps_const, NULL);
            double error = fabs(T_const - T_adaptive);
            fprintf(fp, "%.2f,constant,%.4f,%.8f,%d,%.8f\n", theta0, h, T_const, steps_const, error);
        }
    }

    fclose(fp);
    printf("Arquivo 'analise_completa.csv' gerado com sucesso.\n");
}


/**
Qual o tempo de execu¸c˜ao da simula¸c˜ao para 10 per´ıodos considerando as diferentes es-
trat´egias de passo listadas no item do quadro comparativo?
**/
void run_10_period_analysis() {
    double theta0 = 1.0; // Exemplo com ângulo grande
    double h_vals[] = {0.01, 0.001, 0.0001};
    int n_h = sizeof(h_vals) / sizeof(h_vals[0]);
    double tol = 1e-6; // Tolerância mais alta para o adaptativo
    double h0 = 0.01;
    int num_periods = 10;

    printf("--- Análise de Tempo para %d Períodos (theta0 = %.2f) ---\n", num_periods, theta0);
    printf("method,h,steps,time_s\n");

    // Passos constantes
    for (int j = 0; j < n_h; ++j) {
        double h = h_vals[j];
        int steps_const;
        clock_t tc0 = clock();
        detect_period_constant(theta0, h, num_periods, &steps_const, NULL);
        clock_t tc1 = clock();
        double time_const = (double)(tc1 - tc0) / CLOCKS_PER_SEC;
        printf("constant,%.4f,%d,%.6f\n", h, steps_const, time_const);
    }

    // Integrador adaptativo
    int steps_adapt;
    double T_adapt;
    clock_t ta0 = clock();
    detect_period_adaptive(theta0, tol, h0, num_periods, &T_adapt, &steps_adapt, NULL);
    clock_t ta1 = clock();
    double time_adapt = (double)(ta1 - ta0) / CLOCKS_PER_SEC;
    printf("adaptive,%.4f,%d,%.6f\n", h0, steps_adapt, time_adapt);
}

/**
Baseado no seu experimento, qual o ˆangulo inicial θ0 m´aximo para que a f´ormula simpli-
ficada reporte um per´ıodo com erro menor que 0.001?
**/
void find_max_angle_for_error_threshold() {
    double T_ana = analytic_period();
    double tol = 1e-8; // Usar tolerância bem baixa para ter um T_num preciso
    double h0 = 0.01;
    double error_threshold = 0.001;

    printf("T_analitico = %.8f\n", T_ana);
    
    for (double theta0 = 0.1; theta0 > 0.01; theta0 -= 0.005) {
        double T_num;
        int steps;
        detect_period_adaptive(theta0, tol, h0, 1, &T_num, &steps, NULL);
        double error = fabs(T_num - T_ana);

        printf("theta0 = %.4f, T_num = %.8f, Erro = %.8f\n", theta0, T_num, error);

        if (error < error_threshold) {
            printf("--> Encontrado! Ângulo máximo aproximado para erro < %.4f é %.4f rad.\n", error_threshold, theta0);
            return;
        }
    }
}

/**
Para diferentes valores de θ0, plote o gr´afico de θ × t de um ciclo completo, considerando a
solu¸c˜ao num´erica e a solu¸c˜ao anal´ıtica aproximada.
**/
void generate_plot_data() {
    double theta0_vals[] = {0.1, 0.5, 1.0, 2.0, 3.0}; // Pequeno, médio, grande
    int n_theta = sizeof(theta0_vals) / sizeof(theta0_vals[0]);
    double h = 0.001; // Um h pequeno para o gráfico de passo constante
    double tol = 1e-6;
    int steps;
    double T_num;

    for (int i = 0; i < n_theta; ++i) {
        char filename[100];
        
        // Gerar dados com passo constante
        sprintf(filename, "output/plot_const_theta_%.1f.csv", theta0_vals[i]);
        FILE* fp_const = fopen(filename, "w");
        if (fp_const) {
            fprintf(fp_const, "t,theta\n");
            detect_period_constant(theta0_vals[i], h, 1, &steps, fp_const);
            fclose(fp_const);
        }

        // Gerar dados com passo adaptativo
        sprintf(filename, "output/plot_adapt_theta_%.1f.csv", theta0_vals[i]);
        FILE* fp_adapt = fopen(filename, "w");
        if (fp_adapt) {
            fprintf(fp_adapt, "t,theta\n");
            detect_period_adaptive(theta0_vals[i], tol, h, 1, &T_num, &steps, fp_adapt);
            fclose(fp_adapt);
        }
    }
}