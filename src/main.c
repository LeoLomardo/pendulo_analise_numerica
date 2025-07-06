#include <stdio.h>
#include <math.h>

#include <stdbool.h>
#include "rk.h"
#include "pendulo.h"

int main() {
    double theta0_vals[] = { 0.1, 0.5, 1.0 };
    double h_vals[] = { 0.01, 0.001, 0.0001 };
    int n_theta = sizeof(theta0_vals) / sizeof(theta0_vals[0]);
    int n_h = sizeof(h_vals) / sizeof(h_vals[0]);
    double tol = 1e-6;
    double h0 = 0.01;
    double T_ana = analytic_period();

    // Cabeçalho CSV
    printf("method,theta0,h,T_num,T_ana,steps,time_s,error\n");

    for (int i = 0; i < n_theta; ++i) {
        double theta0 = theta0_vals[i];

        // Passos constantes
        for (int j = 0; j < n_h; ++j) {
            double h = h_vals[j];
            int steps_const;
            clock_t tc0 = clock();
            double T_const = detect_period_constant(theta0, h, &steps_const);
            clock_t tc1 = clock();
            double time_const = (double)(tc1 - tc0) / CLOCKS_PER_SEC;
            double err_const = fabs(T_const - T_ana);
            printf("constant,%.2f,%.4f,%.8f,%.8f,%d,%.6f,%.8f\n",
                   theta0, h, T_const, T_ana, steps_const, time_const, err_const);
        }

        // Integrador adaptativo
        int steps_adapt;
        double T_adapt;
        clock_t ta0 = clock();
        detect_period_adaptive(theta0, tol, h0, &T_adapt, &steps_adapt);
        clock_t ta1 = clock();
        double time_adapt = (double)(ta1 - ta0) / CLOCKS_PER_SEC;
        double err_adapt = fabs(T_adapt - T_ana);
        printf("adaptive,%.2f,%.4f,%.8f,%.8f,%d,%.6f,%.8f\n",
               theta0, h0, T_adapt, T_ana, steps_adapt, time_adapt, err_adapt);
    }

    return 0;
}
