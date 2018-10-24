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
  residuum[0] = x[1];
  residuum[1] = -0.10e0 * x[0] + 0.1e0 * alpha[6] - x[1];
  residuum[2] = x[3];
  residuum[3] = -0.10e0 * x[2] + 0.1e0 * x[0] - x[3];
  residuum[4] = x[5];
  residuum[5] = -0.10e0 * x[4] + 0.1e0 * x[2] - x[5];
  residuum[6] = -w1[1] - omega * w2[0];
  residuum[7] = 0.1e-1 * cos(omega * alpha[0]) * w1[0] + 0.1e-1 * sin(omega * alpha[0]) * w2[0] - 0.1e-1 * cos(omega * (alpha[0] + alpha[3])) * w1[0] - 0.1e-1 * sin(omega * (alpha[0] + alpha[3])) * w2[0] + 0.1e0 * cos(omega * (alpha[0] + alpha[3] + 0.300e3 * sqrt(0.1e1 - exp(-x[0] / 0.90000e5)) + 0.1e1)) * w1[0] + 0.1e0 * sin(omega * (alpha[0] + alpha[3] + 0.300e3 * sqrt(0.1e1 - exp(-x[0] / 0.90000e5)) + 0.1e1)) * w2[0] + w1[1] - omega * w2[1];
  residuum[8] = -w1[3] - omega * w2[2];
  residuum[9] = -cos(omega * alpha[1]) * ((0.1666666667e-1 - 0.8333333333e-3 * alpha[4]) * w1[0] - 0.1e-1 * w1[2]) - sin(omega * alpha[1]) * ((0.1666666667e-1 - 0.8333333333e-3 * alpha[4]) * w2[0] - 0.1e-1 * w2[2]) - 0.1e-1 * cos(omega * (alpha[1] + alpha[4])) * w1[2] - 0.1e-1 * sin(omega * (alpha[1] + alpha[4])) * w2[2] + 0.1e0 * cos(omega * (alpha[1] + alpha[4] + 0.300e3 * sqrt(0.1e1 - exp(-x[2] / 0.90000e5)) + 0.1e1)) * w1[2] + 0.1e0 * sin(omega * (alpha[1] + alpha[4] + 0.300e3 * sqrt(0.1e1 - exp(-x[2] / 0.90000e5)) + 0.1e1)) * w2[2] - cos(omega * (alpha[1] + 0.12e2)) * (0.1e1 / 0.12e2 + 0.8333333333e-3 * alpha[4]) * w1[0] - sin(omega * (alpha[1] + 0.12e2)) * (0.1e1 / 0.12e2 + 0.8333333333e-3 * alpha[4]) * w2[0] + w1[3] - omega * w2[3];
  residuum[10] = -w1[5] - omega * w2[4];
  residuum[11] = -cos(omega * alpha[2]) * ((0.2307692308e-1 - 0.7692307692e-3 * alpha[5]) * w1[2] - 0.1e-1 * w1[4]) - sin(omega * alpha[2]) * ((0.2307692308e-1 - 0.7692307692e-3 * alpha[5]) * w2[2] - 0.1e-1 * w2[4]) - 0.1e-1 * cos(omega * (alpha[2] + alpha[5])) * w1[4] - 0.1e-1 * sin(omega * (alpha[2] + alpha[5])) * w2[4] + 0.1e0 * cos(omega * (alpha[2] + alpha[5] + 0.300e3 * sqrt(0.1e1 - exp(-x[4] / 0.90000e5)) + 0.1e1)) * w1[4] + 0.1e0 * sin(omega * (alpha[2] + alpha[5] + 0.300e3 * sqrt(0.1e1 - exp(-x[4] / 0.90000e5)) + 0.1e1)) * w2[4] - cos(omega * (alpha[2] + 0.13e2)) * (0.1e1 / 0.13e2 + 0.7692307692e-3 * alpha[5]) * w1[2] - sin(omega * (alpha[2] + 0.13e2)) * (0.1e1 / 0.13e2 + 0.7692307692e-3 * alpha[5]) * w2[2] + w1[5] - omega * w2[5];
  residuum[12] = -w2[1] + omega * w1[0];
  residuum[13] = 0.1e-1 * cos(omega * alpha[0]) * w2[0] - 0.1e-1 * sin(omega * alpha[0]) * w1[0] - 0.1e-1 * cos(omega * (alpha[0] + alpha[3])) * w2[0] + 0.1e-1 * sin(omega * (alpha[0] + alpha[3])) * w1[0] + 0.1e0 * cos(omega * (alpha[0] + alpha[3] + 0.300e3 * sqrt(0.1e1 - exp(-x[0] / 0.90000e5)) + 0.1e1)) * w2[0] - 0.1e0 * sin(omega * (alpha[0] + alpha[3] + 0.300e3 * sqrt(0.1e1 - exp(-x[0] / 0.90000e5)) + 0.1e1)) * w1[0] + w2[1] + omega * w1[1];
  residuum[14] = -w2[3] + omega * w1[2];
  residuum[15] = -cos(omega * alpha[1]) * ((0.1666666667e-1 - 0.8333333333e-3 * alpha[4]) * w2[0] - 0.1e-1 * w2[2]) + sin(omega * alpha[1]) * ((0.1666666667e-1 - 0.8333333333e-3 * alpha[4]) * w1[0] - 0.1e-1 * w1[2]) - 0.1e-1 * cos(omega * (alpha[1] + alpha[4])) * w2[2] + 0.1e-1 * sin(omega * (alpha[1] + alpha[4])) * w1[2] + 0.1e0 * cos(omega * (alpha[1] + alpha[4] + 0.300e3 * sqrt(0.1e1 - exp(-x[2] / 0.90000e5)) + 0.1e1)) * w2[2] - 0.1e0 * sin(omega * (alpha[1] + alpha[4] + 0.300e3 * sqrt(0.1e1 - exp(-x[2] / 0.90000e5)) + 0.1e1)) * w1[2] - cos(omega * (alpha[1] + 0.12e2)) * (0.1e1 / 0.12e2 + 0.8333333333e-3 * alpha[4]) * w2[0] + sin(omega * (alpha[1] + 0.12e2)) * (0.1e1 / 0.12e2 + 0.8333333333e-3 * alpha[4]) * w1[0] + w2[3] + omega * w1[3];
  residuum[16] = -w2[5] + omega * w1[4];
  residuum[17] = -cos(omega * alpha[2]) * ((0.2307692308e-1 - 0.7692307692e-3 * alpha[5]) * w2[2] - 0.1e-1 * w2[4]) + sin(omega * alpha[2]) * ((0.2307692308e-1 - 0.7692307692e-3 * alpha[5]) * w1[2] - 0.1e-1 * w1[4]) - 0.1e-1 * cos(omega * (alpha[2] + alpha[5])) * w2[4] + 0.1e-1 * sin(omega * (alpha[2] + alpha[5])) * w1[4] + 0.1e0 * cos(omega * (alpha[2] + alpha[5] + 0.300e3 * sqrt(0.1e1 - exp(-x[4] / 0.90000e5)) + 0.1e1)) * w2[4] - 0.1e0 * sin(omega * (alpha[2] + alpha[5] + 0.300e3 * sqrt(0.1e1 - exp(-x[4] / 0.90000e5)) + 0.1e1)) * w1[4] - cos(omega * (alpha[2] + 0.13e2)) * (0.1e1 / 0.13e2 + 0.7692307692e-3 * alpha[5]) * w2[2] + sin(omega * (alpha[2] + 0.13e2)) * (0.1e1 / 0.13e2 + 0.7692307692e-3 * alpha[5]) * w1[2] + w2[5] + omega * w1[5];
  residuum[18] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) + pow(w2[0], 0.2e1) + pow(w2[1], 0.2e1) + pow(w2[2], 0.2e1) + pow(w2[3], 0.2e1) + pow(w2[4], 0.2e1) + pow(w2[5], 0.2e1) - 0.1e1;
  residuum[19] = w1[0] * w2[0] + w1[1] * w2[1] + w1[2] * w2[2] + w1[3] * w2[3] + w1[4] * w2[4] + w1[5] * w2[5];
}


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* check for proper number of arguments */
 if(nrhs!=6) {
     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nrhs","6 inputs required (some of them are vectors).");
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


 /* create the output matrix */
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)20,mxREAL);

 /* get a pointer to the real data in the output matrix */
 double *residuumPointer = mxGetPr(plhs[0]);

 /* call the computational routine */
 ManifoldEquation(xPointer,alphaPointer,omega,w1Pointer,w2Pointer,residuumPointer);
}