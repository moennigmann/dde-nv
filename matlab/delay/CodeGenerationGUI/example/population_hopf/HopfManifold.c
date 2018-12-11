#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double omega,
  double *w1,
  double *w2,
  double residuum[8])
{
  residuum[0] = alpha[0] * x[1] - alpha[1] * x[0] - alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1];
  residuum[1] = alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - pow(x[1], 0.2e1);
  residuum[2] = cos(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w1[1] + sin(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w2[1] + alpha[1] * w1[0] - (alpha[0] + 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w1[1] - omega * w2[0];
  residuum[3] = -cos(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w1[1] - sin(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w2[1] - (-0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.2e1 * x[1]) * w1[1] - omega * w2[1];
  residuum[4] = cos(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w2[1] - sin(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w1[1] + alpha[1] * w2[0] - (alpha[0] + 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w2[1] + omega * w1[0];
  residuum[5] = -cos(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w2[1] + sin(omega * (0.15e1 - 0.5e0 * exp(-x[1]))) * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w1[1] - (-0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.2e1 * x[1]) * w2[1] + omega * w1[1];
  residuum[6] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w2[0], 0.2e1) + pow(w2[1], 0.2e1) - 0.1e1;
  residuum[7] = w1[0] * w2[0] + w1[1] * w2[1];
}
