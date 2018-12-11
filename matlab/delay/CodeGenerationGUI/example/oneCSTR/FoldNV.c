#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double *w1,
  double *r,
  double residuum[8])
{
  residuum[0] = 0.5e1 / 0.18e2 / alpha[5] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[0]) * w1[0] + 0.5e1 / 0.18e2 / alpha[5] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[1]) * w1[1] + 0.5e1 / 0.18e2 / alpha[5] * (x[5] - x[2]) * w1[2] + (0.5e1 / 0.18e2 / alpha[7] * (x[0] - x[3]) - 0.5e1 / 0.18e2 / alpha[7] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[3])) * w1[3] + (0.5e1 / 0.18e2 / alpha[7] * (x[1] - x[4]) - 0.5e1 / 0.18e2 / alpha[7] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[4])) * w1[4] + 0.5e1 / 0.18e2 / alpha[7] * (x[2] - x[5]) * w1[5] - r[0];
  residuum[1] = -pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) * w1[2] / 0.10e2 - pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) * (x[5] - alpha[6]) * w1[5] / 0.10e2 - r[1];
  residuum[2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * w1[2] - r[2];
  residuum[3] = 0.5e1 / 0.18e2 / alpha[5] * (0.1e1 - x[0]) * w1[0] - 0.5e1 / 0.18e2 / alpha[5] * x[1] * w1[1] + 0.5e1 / 0.18e2 / alpha[5] * (alpha[2] - x[2]) * w1[2] + 0.5e1 / 0.18e2 / alpha[7] * (x[0] - x[3]) * w1[3] + 0.5e1 / 0.18e2 / alpha[7] * (x[1] - x[4]) * w1[4] + 0.5e1 / 0.18e2 / alpha[7] * (x[2] - x[5]) * w1[5] - r[3];
  residuum[4] = alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * w1[2] / 0.10e2 - r[4];
  residuum[5] = (-0.5e1 / 0.18e2 * alpha[3] * pow(alpha[5], -0.2e1) * (0.1e1 - x[0]) - 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[5], -0.2e1) * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[0])) * w1[0] + (0.5e1 / 0.18e2 * alpha[3] * pow(alpha[5], -0.2e1) * x[1] - 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[5], -0.2e1) * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[1])) * w1[1] + (-0.5e1 / 0.18e2 * alpha[3] * pow(alpha[5], -0.2e1) * (alpha[2] - x[2]) - 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[5], -0.2e1) * (x[5] - x[2]) + alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.4e1 / 0.3e1) * (x[2] - alpha[4]) / 0.30e2) * w1[2] - r[5];
  residuum[6] = alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) * w1[5] / 0.10e2 - r[6];
  residuum[7] = (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[7], -0.2e1) * (x[0] - x[3]) + 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[7], -0.2e1) * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[3])) * w1[3] + (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[7], -0.2e1) * (x[1] - x[4]) + 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[7], -0.2e1) * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[4])) * w1[4] + (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[7], -0.2e1) * (x[2] - x[5]) + alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.4e1 / 0.3e1) * (x[5] - alpha[6]) / 0.30e2) * w1[5] - r[7];
}
