#ifndef ODE_H
#define ODE_H

#include <stdio.h>
#define G 9.81
#define L 1.0

int RungeKutta_system_adaptive_h(
    double t0, double t_final, const double y0[], int n_eq,
    void (*f)(double, double[], double[]),
    double tol, double h_initial,
    FILE *outfile, double y_final_out[]);

void rk4_single_step_system(double t, const double y_in[], double h, int n_eq,
                        void (*f)(double, double[], double[]),
                        double y_out[]);

int rk_adaptive_one_step(
    double *t_current, double y_current[], double *h_current, int n_eq,
    void (*f)(double, double[], double[]),
    double tol, double h_min, double h_max);
#endif

