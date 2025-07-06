#ifndef PEN_H
#define PEN_H
#include <time.h>
#include <math.h>
#include <stdbool.h>
// Derivadas para o pêndulo simples: y[0]=θ, y[1]=ω
typedef void (*DerivFunc)(double, double[], double[]);
void f_pendulo(double t, double y[], double dydt[]);
double analytic_period();
double detect_period_constant(double theta0, double h, int *steps_out) ;
int detect_period_adaptive(double theta0, double tol, double h_initial,
                           double *T_num_out, int *steps_out);

#define N_EQ 2

#endif
