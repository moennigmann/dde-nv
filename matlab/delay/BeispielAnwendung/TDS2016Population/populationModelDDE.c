#include "mex.h"
#include <math.h>

void DDErightHandSide (
  double *x,
  double *xtau,
  double *alpha,
  double *p,
  double xdot[2])
{
  xdot[0] = alpha[0] * x[1] - alpha[1] * x[0] - alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * xtau[1];
  xdot[1] = alpha[0] * exp(-alpha[1] * (0.15e1 - 0.5e0 * exp(-x[1]))) * xtau[1] - pow(x[1], 0.2e1);
}

/* The gateway function */                                                                                                                
void mexFunction (int nlhs, mxArray *plhs[],                                                                                              
int nrhs, const mxArray *prhs[])                                                                                                          
{                                                                                                                                         
if(nrhs!=4) {                                                                                                                             
mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nrhs","4 inputs required (state vector, matrix with delayed states and vector with parameters).");
}                                                                                                                                         
if(nlhs==0) {                                                                                                                             
mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nlhs","Please define an output!");                                                                
}                                                                                                                                         
if(nlhs!=1) {                                                                                                                             
mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nlhs","One output required.");                                                                    
}                                                                                                                                         
double *x = mxGetPr(prhs[0]);                                                                                                             
double *xtau = mxGetPr(prhs[1]);                                                                                                          
double *alpha = mxGetPr(prhs[2]);                                                                   
double *p = mxGetPr(prhs[3]);  
/* create the output matrix */                                                                                                            
plhs[0] = mxCreateDoubleMatrix(1,(mwSize)2,mxREAL);                                                                                       
/* get a pointer to the real data in the output matrix */                                                                                 
double *xdot = mxGetPr(plhs[0]);                                                                                                          
DDErightHandSide(x, xtau, alpha, p, xdot);                                                                                                   
}                                                                                                                                         

