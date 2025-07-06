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

// Detecta o período numérico usando passo constante h
double detect_period_constant(double theta0, double h, int *steps_out) {
    double t = 0.0;
    double y[2] = { theta0, 0.0 };
    double y_next[2];
    double prev_omega = y[1];
    double prev_t = t;
    int steps = 0;
    bool first = true;

    while (1) {
        rk4_single_step_system(t, y, h, N_EQ, f_pendulo, y_next);
        steps++;
        double curr_t = t + h;
        double curr_omega = y_next[1];

        if (!first && prev_omega * curr_omega <= 0) {
            double ratio = fabs(prev_omega) / (fabs(prev_omega) + fabs(curr_omega));
            double half_period = prev_t + ratio * h;
            *steps_out = steps;
            return 2.0 * half_period;
        }

        first = false;
        prev_omega = curr_omega;
        prev_t = curr_t;
        t = curr_t;
        y[0] = y_next[0];
        y[1] = y_next[1];
    }
}

// Detecta o período usando integrador adaptativo
int detect_period_adaptive(double theta0, double tol, double h_initial,
                           double *T_num_out, int *steps_out) {
    double t = 0.0;
    double y[2] = { theta0, 0.0 };
    double h = h_initial;
    double h_min = tol * 0.1;
    double h_max = analytic_period() / 4.0;
    double prev_omega = y[1];
    double prev_t = t;
    int steps = 0;
    bool first = true;

    while (1) {
        int status = rk_adaptive_one_step(&t, y, &h, N_EQ, f_pendulo,
                                         tol, h_min, h_max);
        if (status == 0) {
            // Passo rejeitado, tentar de novo sem contar
            continue;
        }
        steps++;
        double curr_t = t;
        double curr_omega = y[1];

        if (!first && prev_omega * curr_omega <= 0) {
            double ratio = fabs(prev_omega) / (fabs(prev_omega) + fabs(curr_omega));
            double half_period = prev_t + ratio * (curr_t - prev_t);
            *T_num_out = 2.0 * half_period;
            *steps_out = steps;
            return 1;
        }

        first = false;
        prev_omega = curr_omega;
        prev_t = curr_t;
    }
}
