#include "mex.h"
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


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

 /* check for proper number of arguments */
 if(nrhs!=4) {
     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldMani:nrhs","4 inputs required (some of them are vectors).");
 }
 if(nlhs==0) {
     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldMani:nlhs","Please define an output!");
 }
 if(nlhs!=1) {
     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldMani:nlhs","One output required.");
 }

 /* get the values of the inputs 
  * scalars are fetched using mxGetScalar and assigned to a type double,
  * vectors are fetched using mxGetPr and assigned to type double pointer
  */

 double *xPointer = mxGetPr(prhs[0]);
 double *alphaPointer = mxGetPr(prhs[1]);
 double *pPointer = mxGetPr(prhs[2]);
 double *wPointer = mxGetPr(prhs[3]);


 /* create the output matrix */
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)
5
,mxREAL);

 /* get a pointer to the real data in the output matrix */
 double *residuumPointer = mxGetPr(plhs[0]);

 /* call the computational routine */
 ManifoldEquation(xPointer, alphaPointer, wPointer, residuumPointer);
}