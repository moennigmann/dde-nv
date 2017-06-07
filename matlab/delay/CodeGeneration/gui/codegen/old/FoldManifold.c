#include <math.h>

void ManifoldEquation (
  double *x,
  double alpha,
  double *w1,
  double residuum[13])
{
  residuum[0] = 0.5e1 / 0.18e2 * alpha4 / alpha6 * (0.1e1 - x[0]) + 0.5e1 / 0.18e2 * alpha1 / alpha6 * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[0]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0];
  residuum[1] = -0.5e1 / 0.18e2 * alpha4 / alpha6 * x[1] + 0.5e1 / 0.18e2 * alpha1 / alpha6 * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[1]) + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0] - 0.2500e4 * exp(-0.7216742844e4 / x[2]) * x[1];
  residuum[2] = 0.5e1 / 0.18e2 * alpha4 / alpha6 * (alpha3 - x[2]) + 0.5e1 / 0.18e2 * alpha1 / alpha6 * (x[5] - x[2]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * x[0] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * x[1] - alpha2 * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha6, -0.1e1 / 0.3e1) * (x[2] - alpha5) / 0.10e2;
  residuum[3] = 0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 * (x[0] - x[3]) - 0.5e1 / 0.18e2 * alpha1 / alpha8 * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[3]);
  residuum[4] = 0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 * (x[1] - x[4]) - 0.5e1 / 0.18e2 * alpha1 / alpha8 * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[4]);
  residuum[5] = 0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 * (x[2] - x[5]) - alpha2 * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha8, -0.1e1 / 0.3e1) * (x[5] - alpha7) / 0.10e2;
  residuum[6] = (-0.5e1 / 0.18e2 * alpha4 / alpha6 - 0.5e1 / 0.18e2 * alpha1 / alpha6 - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w1[0] + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w1[1] + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w1[2] + 0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 * w1[3];
  residuum[7] = (-0.5e1 / 0.18e2 * alpha4 / alpha6 - 0.5e1 / 0.18e2 * alpha1 / alpha6 - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w1[1] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w1[2] + 0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 * w1[4];
  residuum[8] = -0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w1[0] + (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w1[1] + (-0.5e1 / 0.18e2 * alpha4 / alpha6 - 0.5e1 / 0.18e2 * alpha1 / alpha6 + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha2 * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha6, -0.1e1 / 0.3e1) / 0.10e2) * w1[2] + 0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 * w1[5];
  residuum[9] = (-0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 - 0.5e1 / 0.18e2 * alpha1 / alpha8 * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[3] + 0.8333333333e0 * alpha1 / alpha8 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4] + 0.5e1 / 0.18e2 * alpha1 / alpha6 * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[0] - 0.8333333333e0 * alpha1 / alpha6 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[1];
  residuum[10] = 0.4861111111e0 * alpha1 / alpha8 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + (-0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 - 0.5e1 / 0.18e2 * alpha1 / alpha8 * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[4] - 0.4861111111e0 * alpha1 / alpha6 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[0] + 0.5e1 / 0.18e2 * alpha1 / alpha6 * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[1];
  residuum[11] = (-0.1000e4 * (alpha4 / 0.3600e4 + alpha1 / 0.3600e4) / alpha8 - alpha2 * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha8, -0.1e1 / 0.3e1) / 0.10e2) * w1[5] + 0.5e1 / 0.18e2 * alpha1 / alpha6 * w1[2];
  residuum[12] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) - 0.1e1;
}
