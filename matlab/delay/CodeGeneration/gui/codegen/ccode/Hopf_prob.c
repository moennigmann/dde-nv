#include "mex.h"
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double omega,
  double *w1,
  double *w2,
  double residuum[20])
{
  residuum[0] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (0.1e1 - x[0]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[0]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0];
  residuum[1] = -0.5e1 / 0.18e2 * alpha[3] / alpha[5] * x[1] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[1]) + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0] - 0.2500e4 * exp(-0.7216742844e4 / x[2]) * x[1];
  residuum[2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (alpha[2] - x[2]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[5] - x[2]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * x[0] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) / 0.10e2;
  residuum[3] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[0] - x[3]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[3]);
  residuum[4] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[1] - x[4]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[4]);
  residuum[5] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[2] - x[5]) - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) * (x[5] - alpha[6]) / 0.10e2;
  residuum[6] = -cos(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4]) - sin(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4]) - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w1[0] + 0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w1[2] - omega * w2[0];
  residuum[7] = -cos(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[4]) - sin(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[4]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w1[0] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w1[1] - (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w1[2] - omega * w2[1];
  residuum[8] = -0.5e1 / 0.18e2 * cos(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w1[5] - 0.5e1 / 0.18e2 * sin(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w2[5] - 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w1[0] - 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w1[1] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w1[2] - omega * w2[2];
  residuum[9] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[0] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[0] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[7] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4] - omega * w2[3];
  residuum[10] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[1] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[1] - 0.8333333333e0 * alpha[0] / alpha[7] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[4] - omega * w2[4];
  residuum[11] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[2] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[2] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) / 0.10e2) * w1[5] - omega * w2[5];
  residuum[12] = -cos(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4]) + sin(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4]) - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w2[0] + 0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w2[2] + omega * w1[0];
  residuum[13] = -cos(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[4]) + sin(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[4]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w2[0] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w2[1] - (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w2[2] + omega * w1[1];
  residuum[14] = -0.5e1 / 0.18e2 * cos(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w2[5] + 0.5e1 / 0.18e2 * sin(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w1[5] - 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w2[0] - 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w2[1] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w2[2] + omega * w1[2];
  residuum[15] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[0] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[0] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[7] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4] + omega * w1[3];
  residuum[16] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[1] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[1] - 0.8333333333e0 * alpha[0] / alpha[7] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w2[4] + omega * w1[4];
  residuum[17] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[2] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[2] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) / 0.10e2) * w2[5] + omega * w1[5];
  residuum[18] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) + pow(w2[0], 0.2e1) + pow(w2[1], 0.2e1) + pow(w2[2], 0.2e1) + pow(w2[3], 0.2e1) + pow(w2[4], 0.2e1) + pow(w2[5], 0.2e1) - 0.1e1;
  residuum[19] = w1[0] * w2[0] + w1[1] * w2[1] + w1[2] * w2[2] + w1[3] * w2[3] + w1[4] * w2[4] + w1[5] * w2[5];
}
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double omega,
  double *w1,
  double *w2,
  double residuum[20])
{
  residuum[0] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (0.1e1 - x[0]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[0]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0];
  residuum[1] = -0.5e1 / 0.18e2 * alpha[3] / alpha[5] * x[1] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[1]) + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0] - 0.2500e4 * exp(-0.7216742844e4 / x[2]) * x[1];
  residuum[2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (alpha[2] - x[2]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[5] - x[2]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * x[0] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) / 0.10e2;
  residuum[3] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[0] - x[3]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[3]);
  residuum[4] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[1] - x[4]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[4]);
  residuum[5] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[2] - x[5]) - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) * (x[5] - alpha[6]) / 0.10e2;
  residuum[6] = -cos(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4]) - sin(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4]) - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w1[0] + 0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w1[2] - omega * w2[0];
  residuum[7] = -cos(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[4]) - sin(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[4]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w1[0] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w1[1] - (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w1[2] - omega * w2[1];
  residuum[8] = -0.5e1 / 0.18e2 * cos(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w1[5] - 0.5e1 / 0.18e2 * sin(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w2[5] - 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w1[0] - 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w1[1] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w1[2] - omega * w2[2];
  residuum[9] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[0] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[0] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[7] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4] - omega * w2[3];
  residuum[10] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[1] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[1] - 0.8333333333e0 * alpha[0] / alpha[7] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[4] - omega * w2[4];
  residuum[11] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[2] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[2] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) / 0.10e2) * w1[5] - omega * w2[5];
  residuum[12] = -cos(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4]) + sin(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4]) - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w2[0] + 0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w2[2] + omega * w1[0];
  residuum[13] = -cos(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[4]) + sin(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[4]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w2[0] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w2[1] - (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w2[2] + omega * w1[1];
  residuum[14] = -0.5e1 / 0.18e2 * cos(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w2[5] + 0.5e1 / 0.18e2 * sin(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w1[5] - 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w2[0] - 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w2[1] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w2[2] + omega * w1[2];
  residuum[15] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[0] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[0] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[7] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4] + omega * w1[3];
  residuum[16] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[1] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[1] - 0.8333333333e0 * alpha[0] / alpha[7] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w2[4] + omega * w1[4];
  residuum[17] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[2] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[2] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) / 0.10e2) * w2[5] + omega * w1[5];
  residuum[18] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) + pow(w2[0], 0.2e1) + pow(w2[1], 0.2e1) + pow(w2[2], 0.2e1) + pow(w2[3], 0.2e1) + pow(w2[4], 0.2e1) + pow(w2[5], 0.2e1) - 0.1e1;
  residuum[19] = w1[0] * w2[0] + w1[1] * w2[1] + w1[2] * w2[2] + w1[3] * w2[3] + w1[4] * w2[4] + w1[5] * w2[5];
}
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double omega,
  double *w1,
  double *w2,
  double residuum[20])
{
  residuum[0] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (0.1e1 - x[0]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[0]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0];
  residuum[1] = -0.5e1 / 0.18e2 * alpha[3] / alpha[5] * x[1] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[1]) + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0] - 0.2500e4 * exp(-0.7216742844e4 / x[2]) * x[1];
  residuum[2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (alpha[2] - x[2]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[5] - x[2]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * x[0] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) / 0.10e2;
  residuum[3] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[0] - x[3]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 * x[3] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[3]);
  residuum[4] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[1] - x[4]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (x[4] / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - x[4]);
  residuum[5] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * (x[2] - x[5]) - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) * (x[5] - alpha[6]) / 0.10e2;
  residuum[6] = -cos(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4]) - sin(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4]) - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w1[0] + 0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w1[2] - omega * w2[0];
  residuum[7] = -cos(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[4]) - sin(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[4]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w1[0] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w1[1] - (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w1[2] - omega * w2[1];
  residuum[8] = -0.5e1 / 0.18e2 * cos(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w1[5] - 0.5e1 / 0.18e2 * sin(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w2[5] - 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w1[0] - 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w1[1] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w1[2] - omega * w2[2];
  residuum[9] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[0] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[0] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[7] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4] - omega * w2[3];
  residuum[10] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[1] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[1] - 0.8333333333e0 * alpha[0] / alpha[7] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w1[4] - omega * w2[4];
  residuum[11] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[2] - 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[2] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) / 0.10e2) * w1[5] - omega * w2[5];
  residuum[12] = -cos(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4]) + sin(0.1750e4 * omega / alpha[0]) * (0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[3] - 0.4861111111e0 * alpha[0] / alpha[5] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[4]) - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w2[0] + 0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w2[2] + omega * w1[0];
  residuum[13] = -cos(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w2[4]) + sin(0.1750e4 * omega / alpha[0]) * (-0.8333333333e0 * alpha[0] / alpha[5] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w1[3] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1)) * w1[4]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w2[0] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w2[1] - (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w2[2] + omega * w1[1];
  residuum[14] = -0.5e1 / 0.18e2 * cos(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w2[5] + 0.5e1 / 0.18e2 * sin(0.1750e4 * omega / alpha[0]) * alpha[0] / alpha[5] * w1[5] - 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w2[0] - 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w2[1] - (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w2[2] + omega * w1[2];
  residuum[15] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[0] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[0] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.35e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.1050e2 * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w2[3] - 0.4861111111e0 * alpha[0] / alpha[7] * x[3] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[4] + omega * w1[3];
  residuum[16] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[1] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[1] - 0.8333333333e0 * alpha[0] / alpha[7] * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) * w2[3] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - 0.5e1 / 0.18e2 * alpha[0] / alpha[7] * (0.1e1 / (0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0) - 0.5e0 * x[4] * pow(0.30e1 * x[3] + 0.5e0 * x[4] + 0.5e0, -0.2e1) - 0.1e1)) * w2[4] + omega * w1[4];
  residuum[17] = -0.1000e4 * cos(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w2[2] + 0.1000e4 * sin(0.1750e4 * omega / (alpha[0] + alpha[3])) * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] * w1[2] - (-0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[7], -0.1e1 / 0.3e1) / 0.10e2) * w2[5] + omega * w1[5];
  residuum[18] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) + pow(w2[0], 0.2e1) + pow(w2[1], 0.2e1) + pow(w2[2], 0.2e1) + pow(w2[3], 0.2e1) + pow(w2[4], 0.2e1) + pow(w2[5], 0.2e1) - 0.1e1;
  residuum[19] = w1[0] * w2[0] + w1[1] * w2[1] + w1[2] * w2[2] + w1[3] * w2[3] + w1[4] * w2[4] + w1[5] * w2[5];
}
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

/* check for proper number of arguments */
 if(nrhs!=12) {
     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nrhs","12 inputs required (some of them are vectors).");
}
 if(nlhs==0) {
     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nlhs","Please define an output!");
}
 if(nlhs!=1) {
     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nlhs","One output required.");
}

 /* get the values of the inputs
  * scalars are fetched using mxGetScalar and assigned to a type double,
  * vectors are fetched using mxGetPr and assigned to type double pointer
  */

 double *xPointer = mxGetPr(prhs[0]);
 double *alphaPointer = mxGetPr(prhs[1]);
 double *pPointer = mxGetPr(prhs[2]);
 double omega = mxGetScalar(prhs[3]);
 double *w1Pointer = mxGetPr(prhs[4]);
 double *w2Pointer = mxGetPr(prhs[5]);
 double *v1Pointer = mxGetPr(prhs[6]);
 double *v2Pointer = mxGetPr(prhs[7]);
 double g1 = mxGetScalar(prhs[8]);
 double g2 = mxGetScalar(prhs[9]);
 double *uPointer = mxGetPr(prhs[10]);
 double *rPointer = mxGetPr(prhs[11]);



 /* create the output matrix */
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)47,mxREAL);

 /* get a pointer to the real data in the output matrix */
 double *residuumPointer = mxGetPr(plhs[0]);

 /* call the computational routine */
 ManifoldEquation(xPointer,alphaPointer,omega,w1Pointer,w2Pointer,v1Pointer,
         v2Pointer,g1,g2,uPointer,rPointer,residuumPointer);
}
