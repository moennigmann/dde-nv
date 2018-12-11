#include "mex.h"
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
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)8,mxREAL);

 /* get a pointer to the real data in the output matrix */
 double *residuumPointer = mxGetPr(plhs[0]);

 /* call the computational routine */
 ManifoldEquation(xPointer,alphaPointer,omega,w1Pointer,w2Pointer,residuumPointer);
}
