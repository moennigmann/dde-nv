#include "mex.h"
#include <math.h>

void DDErightHandSide (
  double *x,
  double *xtau,
  double *alpha,
  double *p,
  double xdot[6])
{
  xdot[0] = 0.1000e4 * p[0] * x[1] + 0.1000e4 * x[2] * x[0] - 0.5000e4 * x[2] * x[1] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[0] + 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[1] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[3] + 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[4];
  xdot[1] = -0.1000e4 * x[0] * p[0] + 0.5000e4 * x[2] * x[0] + 0.1000e4 * x[2] * x[1] - 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[0] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[1] - 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[3] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[4];
  xdot[2] = 0.1000e1 * alpha[0] - 0.1000e1 * x[2] - 0.1000e1 * (0.2e1 * x[2] + 0.1e1) * (pow(x[0], 0.2e1) + pow(x[1], 0.2e1));
  xdot[3] = 0.1000e4 * p[0] * x[4] + 0.1000e4 * x[5] * x[3] - 0.5000e4 * x[5] * x[4] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[3] + 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[4] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[0] + 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[1];
  xdot[4] = -0.1000e4 * p[0] * x[3] + 0.5000e4 * x[5] * x[3] + 0.1000e4 * x[5] * x[4] - 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[3] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[4] - 0.2500000000e0 * sin(0.1000e4 * p[0] + 0.1e1) * xtau[0] + 0.2500000000e0 * cos(0.1000e4 * p[0] + 0.1e1) * xtau[1];
  xdot[5] = 0.1000e1 * alpha[1] - 0.1000e1 * x[5] - 0.1000e1 * (0.2e1 * x[5] + 0.1e1) * (pow(x[3], 0.2e1) + pow(x[4], 0.2e1));
}


/* The gateway function */                                                                                                                
void mexFunction (int nlhs, mxArray *plhs[],                                                                                              
int nrhs, const mxArray *prhs[])                                                                                                          
{                                                                                                                                         
if(nrhs!=4) {                                                                                                                             
mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nrhs","4 inputs required (state vector, matrix with delayed states and vector with parameters (uncertain und certain)).");
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
plhs[0] = mxCreateDoubleMatrix(1,(mwSize)6,mxREAL);                                                                                       
/* get a pointer to the real data in the output matrix */                                                                                 
double *xdot = mxGetPr(plhs[0]);                                                                                                          
DDErightHandSide(x, xtau, alpha, p, xdot);                                                                                                   
}                                                                                                                                         

