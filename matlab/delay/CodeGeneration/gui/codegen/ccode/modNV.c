#include "mex.h"
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


/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{

 /* check for proper number of arguments */
 if(nrhs!=8) {  
     mexErrMsgIdAndTxt("MyToolbox:populationModelModFoldNVPop:nrhs","8 inputs required (some of them are vectors).");
 }
 if(nlhs==0) { 
    mexErrMsgIdAndTxt("MyToolbox:populationModelModFoldNVPop:nlhs","Please define an output!");
 }
 if(nlhs!=1) {
     mexErrMsgIdAndTxt("MyToolbox:populationModelModFoldNVPop:nlhs","One output required.");
 }

 /* get the values of the inputs 
  * scalars are fetched using mxGetScalar and assigned to a type double,
  * vectors are fetched using mxGetPr and assigned to type double pointer
  */

 double *xPointer = mxGetPr(prhs[0]);
 double *alphaPointer = mxGetPr(prhs[1]);
 double *pPointer = mxGetPr(prhs[2]);
 double *wPointer = mxGetPr(prhs[3]);
 double *vPointer = mxGetPr(prhs[4]);
 double g1 = mxGetScalar(prhs[5]);
 double *uPointer = mxGetPr(prhs[6]);
 double *rPointer = mxGetPr(prhs[7]);


 /* create the output matrix */
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)7,mxREAL);

 /* get a pointer to the real data in the output matrix */
 double *residuumPointer = mxGetPr(plhs[0]);

 /* call the computational routine */
 ManifoldEquation(xPointer,alphaPointer, pPointer, wPointer,vPointer,g1,uPointer,rPointer, residuumPointer);
}