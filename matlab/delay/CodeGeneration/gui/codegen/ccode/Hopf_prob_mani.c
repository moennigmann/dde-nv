#include "mex.h"
void ManifoldEquation (
  double x,
  double alpha,
  double w1,
  double residuum[0])
{
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
