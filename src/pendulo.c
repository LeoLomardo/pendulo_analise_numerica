#include "pendulo.h"
#include "rk.h"

void f_pendulo(double t, double y[], double dydt[]) {
    dydt[0] = y[1];
    dydt[1] = -(G/L) * sin(y[0]);
}

// Período analítico (pequeno ângulo)
double analytic_period() {
    return 2.0 * M_PI * sqrt(L / G);
}


// Detecta o período numérico usando passo constante h.
double detect_period_constant(double theta0, double h, int num_periods, int *steps_out, FILE* outfile) {
    double t = 0.0;
    double y[2] = { theta0, 0.0 };
    double y_next[2];
    double prev_omega = y[1];
    int steps = 0;
    
    int zero_crossings = 0;
    double first_half_period_time = 0.0;
    double total_time = 0.0;

    // Salva o ponto inicial.
    if (outfile) {
        fprintf(outfile, "%.8f,%.8f\n", t, y[0]);
    }

    while (zero_crossings < 2 * num_periods) {
        rk4_single_step_system(t, y, h, N_EQ, f_pendulo, y_next);
        steps++;
        double curr_t = t + h;
        double curr_omega = y_next[1];

        // Salva o ponto atual.
        if (outfile) {
            fprintf(outfile, "%.8f,%.8f\n", curr_t, y[0]);
        }

        // Detecta cruzamento por zero da velocidade angular (inversão de movimento)
        if (prev_omega * curr_omega <= 0 && t > 0) {
            zero_crossings++;
            
            // Interpolação linear para encontrar o tempo exato do cruzamento
            double ratio = fabs(prev_omega) / (fabs(prev_omega) + fabs(curr_omega));
            double interpolated_time = t + ratio * h;

            if (zero_crossings == 1) { // Primeiro meio período
                first_half_period_time = interpolated_time;
            }
            // Ao final de `num_periods` (2*num_periods meios-períodos)
            if (zero_crossings == 2 * num_periods) {
                 total_time = interpolated_time;
            }
        }

        prev_omega = curr_omega;
        t = curr_t;
        y[0] = y_next[0];
        y[1] = y_next[1];
    }
    
    *steps_out = steps;
   
    // Retorna o período médio
    return (2.0 * total_time) / (double)(2 * num_periods);
}


//  Detecta o período usando integrador adaptativo.
int detect_period_adaptive(double theta0, double tol, double h_initial, int num_periods,
                           double *T_num_out, int *steps_out, FILE* outfile) {
    double t = 0.0;
    double y[2] = { theta0, 0.0 };
    double h = h_initial;
    double h_min = tol * 0.1;
    double h_max = analytic_period() / 4.0;
    double prev_omega = y[1];
    double prev_t = t;
    int steps = 0;

    int zero_crossings = 0;
    double total_time = 0.0;
    
    if (outfile) {
        fprintf(outfile, "%.8f,%.8f\n", t, y[0]);
    }

    while (zero_crossings < 2 * num_periods) {
        double t_before_step = t;
        int status = rk_adaptive_one_step(&t, y, &h, N_EQ, f_pendulo,
                                         tol, h_min, h_max);
        if (status == 0) { // Passo rejeitado
            continue;
        }
        steps++;
        
        if (outfile) {
            fprintf(outfile, "%.8f,%.8f\n", t, y[0]);
        }

        double curr_omega = y[1];
        // Detecta cruzamento por zero
        if (prev_omega * curr_omega <= 0 && t_before_step > 0) {
            zero_crossings++;

            double ratio = fabs(prev_omega) / (fabs(prev_omega) + fabs(curr_omega));
            double interpolated_time = prev_t + ratio * (t - prev_t);
            
            if (zero_crossings == 2 * num_periods) {
                total_time = interpolated_time;
            }
        }

        prev_omega = curr_omega;
        prev_t = t;
    }

    *steps_out = steps;
    // Retorna o período médio
    *T_num_out = (2.0 * total_time) / (double)(2 * num_periods);
    return 1;
}