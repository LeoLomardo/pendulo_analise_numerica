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

#endif