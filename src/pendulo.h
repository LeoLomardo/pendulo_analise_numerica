#ifndef PENDULO_H
#define PENDULO_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// Constantes Físicas e do Sistema
#define G 9.81      // Aceleração da gravidade
#define L 1.0       // Comprimento do pêndulo
#define N_EQ 2      // Número de equações (theta, omega)

/**
 * @brief Define o sistema de EDOs para o pêndulo.
 * y[0] = theta, y[1] = omega
 */
void f_pendulo(double t, double y[], double dydt[]);

/**
 * @brief Calcula o período analítico para pequenas oscilações.
 */
double analytic_period();

/**
 * @brief Detecta o período numérico usando passo constante h.
 * @param theta0 Ângulo inicial.
 * @param h Tamanho do passo.
 * @param num_periods Número de períodos a simular.
 * @param steps_out Ponteiro para armazenar o número total de passos.
 * @param outfile Ponteiro de arquivo para salvar dados de theta vs t (pode ser NULL).
 * @return O período médio calculado ao longo de num_periods.
 */
double detect_period_constant(double theta0, double h, int num_periods, int *steps_out, FILE* outfile);

/**
 * @brief Detecta o período numérico usando passo adaptativo.
 * @param theta0 Ângulo inicial.
 * @param tol Tolerância de erro.
 * @param h_initial Tamanho inicial do passo.
 * @param num_periods Número de períodos a simular.
 * @param T_num_out Ponteiro para armazenar o período médio final.
 * @param steps_out Ponteiro para armazenar o número total de passos.
 * @param outfile Ponteiro de arquivo para salvar dados de theta vs t (pode ser NULL).
 * @return 1 em caso de sucesso, 0 em caso de falha.
 */
int detect_period_adaptive(double theta0, double tol, double h_initial, int num_periods, double *T_num_out, int *steps_out, FILE* outfile);


#endif