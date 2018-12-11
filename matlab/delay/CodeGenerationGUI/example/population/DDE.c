#include <math.h>

void DDErightHandSide (
  double *x,
  double *xtau,
  double *alpha,
  double xdot[2])
{
  xdot[0] = alpha[0] * x[1] - alpha[1] * x[0] - alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * xtau[1];
  xdot[1] = alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * xtau[1] - pow(x[1], 0.2e1);
}
