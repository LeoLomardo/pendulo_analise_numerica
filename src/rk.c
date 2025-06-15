#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Realiza um único passo do método Runge-Kutta de 4ª ordem para um sistema de EDOs.
 * y_out = y_in + resultado_do_passo_rk4
 * @param t Tempo atual.
 * @param y_in Vetor do estado atual [y1, y2, ..., yn].
 * @param h Tamanho do passo.
 * @param n_eq Número de equações no sistema.
 * @param f Ponteiro para a função de derivadas f(t, y, dydt).
 * @param y_out Vetor onde o resultado do passo y(t+h) será armazenado.
 */
void rk4_single_step_system(double t, const double y_in[], double h, int n_eq,
                            void (*f)(double, double[], double[]),
                            double y_out[])
{
    double k1[n_eq], k2[n_eq], k3[n_eq], k4[n_eq];
    double y_temp[n_eq];
    int i;

    // Calcular k1
    f(t, (double *)y_in, k1); // k1 = f(t, y_in)

    // Calcular k2
    for (i = 0; i < n_eq; ++i)
    {
        y_temp[i] = y_in[i] + (h / 2.0) * k1[i];
    }
    f(t + h / 2.0, y_temp, k2); // k2 = f(t + h/2, y_in + h/2 * k1)

    // Calcular k3
    for (i = 0; i < n_eq; ++i)
    {
        y_temp[i] = y_in[i] + (h / 2.0) * k2[i];
    }
    f(t + h / 2.0, y_temp, k3); // k3 = f(t + h/2, y_in + h/2 * k2)

    // Calcular k4
    for (i = 0; i < n_eq; ++i)
    {
        y_temp[i] = y_in[i] + h * k3[i];
    }
    f(t + h, y_temp, k4); // k4 = f(t + h, y_in + h * k3)

    // Combinar para o resultado final
    for (i = 0; i < n_eq; ++i)
    {
        y_out[i] = y_in[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

/**
 * @brief Realiza um passo de integração com controle de erro para RK4 adaptativo.
 *        Baseado na estratégia de dobrar o passo.
 * @param t_current Ponteiro para o tempo atual (será atualizado).
 * @param y_current Vetor do estado atual (será atualizado).
 * @param h_current Ponteiro para o tamanho do passo atual (será atualizado para o próximo passo).
 * @param n_eq Número de equações.
 * @param f Ponteiro para a função de derivadas.
 * @param tol Tolerância de erro desejada (para o erro local em y[0]).
 * @param h_min Limite inferior para o tamanho do passo.
 * @param h_max Limite superior para o tamanho do passo.
 * @return 1 se o passo foi aceito, 0 se foi rejeitado (e h_current foi reduzido).
 */
int rk_adaptive_one_step(
    double *t_current, double y_current[], double *h_current, int n_eq,
    void (*f)(double, double[], double[]),
    double tol, double h_min, double h_max)
{
    double y1[n_eq];       // Resultado de um passo h
    double y2_half1[n_eq]; // Resultado do primeiro meio passo (h/2)
    double y2[n_eq];       // Resultado de dois meio passos (h/2 + h/2)

    double h = *h_current;
    int i;

    // 1. Calcular y1 (um passo de h)
    rk4_single_step_system(*t_current, y_current, h, n_eq, f, y1);

    // 2. Calcular y2 (dois passos de h/2)
    double h_half = h / 2.0;
    rk4_single_step_system(*t_current, y_current, h_half, n_eq, f, y2_half1);
    rk4_single_step_system(*t_current + h_half, y2_half1, h_half, n_eq, f, y2);

    // 3. Estimar o erro (truncamento local)
    // O erro é estimado como |y2[0] - y1[0]| / 15.0 para a componente theta (y[0])
    // (para RK4, o denominador é 2^p - 1 = 2^4 - 1 = 15)
    double error_estimate = fabs(y2[0] - y1[0]) / 15.0; // Foco em theta

    // 4. Decidir se aceita o passo e calcular o novo h
    double safety_factor = 0.9; // Fator de segurança < 1
    double h_new;

    if (error_estimate <= tol || h <= h_min * 1.0001)
    { // Passo aceito (ou h já é mínimo)
        *t_current += h;
        for (i = 0; i < n_eq; ++i)
        {
            // Usar a solução mais precisa (y2) e adicionar a estimativa de erro
            // para obter uma solução de ordem superior (local extrapolation)
            y_current[i] = y2[i] + (y2[i] - y1[i]) / 15.0;
        }

        // Calcular h_new para o próximo passo
        if (error_estimate == 0.0)
        { // Aumenta o passo (ex: dobrar)
            h_new = h * 2.0;
        }
        else
        {
            // safety_factor * (desejado/obtido)^(1/5)
            h_new = h * safety_factor * pow(tol / error_estimate, 0.20); // 0.20 = 1/(p+1) para RK4 se p=4 ou 1/p
                                                                         // A literatura geralmente usa 1/ (ordem do erro + 1) -> 1/5 para RK4
        }
        *h_current = fmin(fmax(h_new, h_min), h_max); // Limitar h_new
        return 1;                                     // Passo aceito
    }
    else
    { // Erro muito grande, rejeitar o passo e reduzir h
        h_new = h * safety_factor * pow(tol / error_estimate, 0.20);
        *h_current = fmax(h_new, h_min); // Reduzir h, mas não abaixo de h_min
        // Não avançar *t_current nem y_current
        return 0; // Passo rejeitado
    }
}

/**
 * @brief Resolve um sistema de EDOs usando Runge-Kutta de 4ª ordem com passo adaptativo.
 * @param t0 Tempo inicial.
 * @param t_final Tempo final.
 * @param y0 Vetor de condições iniciais.
 * @param n_eq Número de equações.
 * @param f Ponteiro para a função de derivadas.
 * @param tol Tolerância de erro desejada.
 * @param h_initial Estimativa inicial para o tamanho do passo.
 * @param outfile Ponteiro para o arquivo de saída (NULL se não quiser salvar).
 * @param y_final_out Vetor para armazenar o estado final (opcional, pode ser NULL).
 * @return Número de passos aceitos.
 */
int RungeKutta_system_adaptive_h(
    double t0, double t_final, const double y0[], int n_eq,
    void (*f)(double, double[], double[]),
    double tol, double h_initial,
    FILE *outfile, double y_final_out[])
{
    double t = t0;
    double y[n_eq];
    double h = h_initial;
    int i;

    double h_min = tol * 0.1;             // Exemplo de h_min, pode ser ajustado
    double h_max = (t_final - t0) / 10.0; // Exemplo de h_max

    int accepted_steps = 0;
    int rejected_steps = 0; // Para estatísticas, se desejar

    for (i = 0; i < n_eq; ++i)
        y[i] = y0[i];

    if (outfile)
    { // Salvar ponto inicial
        fprintf(outfile, "%.8e", t);
        for (i = 0; i < n_eq; ++i)
            fprintf(outfile, " %.8e", y[i]);
        fprintf(outfile, "\n");
    }

    while (t < t_final)
    {
        if (t + h > t_final)
        { // Ajustar o passo para não ultrapassar t_final
            h = t_final - t;
        }
        if (h <= 1e-12)
            break; // Evitar passo efetivamente nulo

        double y_temp_before_step[n_eq]; // Salvar estado caso o passo seja rejeitado
        memcpy(y_temp_before_step, y, n_eq * sizeof(double));
        double t_before_step = t;
        double h_try = h; // h que será tentado (e possivelmente modificado por rk_adaptive_one_step)

        int status = rk_adaptive_one_step(&t, y, &h_try, n_eq, f, tol, h_min, h_max);

        if (status == 1)
        { // Passo aceito
            accepted_steps++;
            h = h_try; // h_try foi atualizado para o próximo passo sugerido
            if (outfile)
            {
                fprintf(outfile, "%.8e", t);
                for (i = 0; i < n_eq; ++i)
                    fprintf(outfile, " %.8e", y[i]);
                fprintf(outfile, "\n");
            }
        }
        else
        { // Passo rejeitado
            rejected_steps++;
            memcpy(y, y_temp_before_step, n_eq * sizeof(double)); // Restaurar y
            t = t_before_step;                                    // Restaurar t
            h = h_try;                                            // h_try foi reduzido pela rk_adaptive_one_step
            // O loop continuará e tentará novamente com o novo h reduzido.
        }
    }
    if (y_final_out)
    {
        for (i = 0; i < n_eq; ++i)
            y_final_out[i] = y[i];
    }
    // printf("Adaptativo: Aceitos=%d, Rejeitados=%d\n", accepted_steps, rejected_steps);
    return accepted_steps;
}