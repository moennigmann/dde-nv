#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double *w,
  double residuum[5])
{
  residuum[0] = alpha[0] * x[1] - alpha[1] * x[0] - alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1];
  residuum[1] = alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - pow(x[1], 0.2e1);
  residuum[2] = exp(0.15e0 - 0.5e-1 * exp(-x[1])) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] + alpha[1] * w[0] - (alpha[0] + 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w[1] - 0.1e0 * w[0];
  residuum[3] = -exp(0.15e0 - 0.5e-1 * exp(-x[1])) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] - (-0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.2e1 * x[1]) * w[1] - 0.1e0 * w[1];
  residuum[4] = pow(w[0], 0.2e1) + pow(w[1], 0.2e1) - 0.1e1;
}
