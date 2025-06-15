#include "rk.h"
#include <stdio.h>
#include <math.h>

static int N;
static void f(double t, double y[], double dydt[])
{
  N++;
  dydt[0] = t * y[0] + t * t * t;
}

#define S 10.0542731796122

int main(void)
{
  double h = 0.001;
  double t1 = 2.4, t0 = 0, y0 = -1;
  double y_final[1];
  double tol = 1e-12;

  printf("----------------------------------------------------\n");

  N = 0;
  int passos = RungeKutta_system_adaptive_h(t0, t1, &y0, 1, f, tol, h, NULL, y_final);
  printf("RungeKuttaAdapt: N=%d rk=%g diff=%g\n", N, y_final[0], fabs(S - y_final[0]));

  printf("----------------------------------------------------\n");
  return 0;
}
