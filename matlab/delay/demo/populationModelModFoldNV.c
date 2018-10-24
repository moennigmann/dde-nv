#include "mex.h"
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double *p,
  double *w,
  double *v,
  double g1,
  double *u,
  double *r,
  double residuum[7])
{
  residuum[0] = alpha[1] * v[0] - 0.1e0 * v[0] + 0.2e1 * g1 * w[0];
  residuum[1] = -exp(0.15e0 - 0.5e-1 * exp(-x[1])) * (-alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * v[0] + alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * v[1]) - (alpha[0] + 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * v[0] - (-0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.2e1 * x[1]) * v[1] - 0.1e0 * v[1] + 0.2e1 * g1 * w[1];
  residuum[2] = -alpha[1] * u[0];
  residuum[3] = -alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * u[0] + alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * u[1] - exp(0.15e0 - 0.5e-1 * exp(-x[1])) * (0.5e0 * v[0] * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] - 0.5e0 * v[1] * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1]) + exp(0.15e0 - 0.5e-1 * exp(-x[1])) * (0.5e-1 * v[0] * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] * exp(-x[1]) - 0.5e-1 * v[1] * alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] * exp(-x[1])) + (alpha[0] + 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * u[0] + (-0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.2e1 * x[1]) * u[1] - v[0] * (-0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.25e0 * alpha[0] * pow(alpha[1], 0.2e1) * pow(exp(-x[1]), 0.2e1) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w[1] - v[1] * (0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] + 0.25e0 * alpha[0] * pow(alpha[1], 0.2e1) * pow(exp(-x[1]), 0.2e1) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.2e1) * w[1];
  residuum[4] = -exp(0.15e0 - 0.5e-1 * exp(-x[1])) * (-v[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] + v[1] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1]) + (x[1] - exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * u[0] + exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] * u[1] - v[0] * (0.1e1 + 0.5e0 * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w[1] + 0.5e0 * v[1] * alpha[1] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] * w[1] - r[0];
  residuum[5] = -exp(0.15e0 - 0.5e-1 * exp(-x[1])) * (-v[0] * alpha[0] * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1] + v[1] * alpha[0] * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * w[1]) + (-x[0] - alpha[0] * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * u[0] + alpha[0] * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] * u[1] - v[0] * (-w[0] + (0.5e0 * alpha[0] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] + 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w[1]) - v[1] * (-0.5e0 * alpha[0] * exp(-x[1]) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1] - 0.5e0 * alpha[0] * alpha[1] * exp(-x[1]) * (-0.15e1 + 0.5e0 * exp(-x[1])) * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * x[1]) * w[1] - r[1];
  residuum[6] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) - 0.1e1;
}

/* The gateway function */                                                                                           
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])                                         
{                                                                                                                    
                                                                                                                     
 /* check for proper number of arguments */                                                                          
 if(nrhs!=9) {                                                                                                       
     mexErrMsgIdAndTxt("MyToolbox:populationModelModFoldNVPop:nrhs","9 inputs required (some of them are vectors).");
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
 double *kPointer = mxGetPr(prhs[6]);                                                                                
 double *uPointer = mxGetPr(prhs[7]);                                                                                
 double *rPointer = mxGetPr(prhs[8]);                                                                                
                                                                                                                     
                                                                                                                     
 /* create the output matrix */                                                                                      
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)7,mxREAL);                                                                 
                                                                                                                     
 /* get a pointer to the real data in the output matrix */                                                           
 double *residuumPointer = mxGetPr(plhs[0]);                                                                         
                                                                                                                     
 /* call the computational routine */                                                                                
 ManifoldEquation(xPointer,alphaPointer, pPointer, wPointer,vPointer,g1,uPointer,rPointer, residuumPointer);                   
}                                                                                                                    

