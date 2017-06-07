#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double w,
  double v,
  double g1,
  double u,
  double *r,
  double residuum[2])
{
  residuum[0] = (x[1] - exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w1[0] + exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] * w1[1] - r[0];
  residuum[1] = (-x[0] - alpha[0] * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w1[0] + alpha[0] * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] * w1[1] - r[1];
}
