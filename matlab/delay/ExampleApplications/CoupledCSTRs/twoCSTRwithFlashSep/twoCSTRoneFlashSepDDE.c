#include "mex.h"
#include <math.h>

void DDErightHandSide (
  double *x,
  double *xtau,
  double *alpha,
  double *p,
  double xdot[9])
{
  xdot[0] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (0.1e1 - x[0]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 * xtau[6] / (0.30e1 * xtau[6] + 0.5e0 * xtau[7] + 0.5e0) - x[0]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0];
  xdot[1] = -0.5e1 / 0.18e2 * alpha[3] / alpha[5] * x[1] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (xtau[7] / (0.30e1 * xtau[6] + 0.5e0 * xtau[7] + 0.5e0) - x[1]) + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0] - 0.2500e4 * exp(-0.7216742844e4 / x[2]) * x[1];
  xdot[2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (alpha[2] - x[2]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (xtau[8] - x[2]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * x[0] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) / 0.100e3;
  xdot[3] = 0.5e1 / 0.18e2 * alpha[7] / alpha[9] * (0.1e1 - x[3]) + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * (xtau[9] - x[3]) - 0.2770e4 * exp(-0.6013952370e4 / x[5]) * x[3];
  xdot[4] = -0.5e1 / 0.18e2 * alpha[7] / alpha[9] * x[4] + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * (xtau[10] - x[4]) + 0.2770e4 * exp(-0.6013952370e4 / x[5]) * x[3] - 0.2500e4 * exp(-0.7216742844e4 / x[5]) * x[4];
  xdot[5] = 0.5e1 / 0.18e2 * alpha[7] / alpha[9] * (alpha[6] - x[5]) + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * (xtau[11] - x[5]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[5]) * x[3] + 0.4166666675e8 * exp(-0.7216742844e4 / x[5]) * x[4] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[9], -0.1e1 / 0.3e1) * (x[5] - alpha[8]) / 0.100e3;
  xdot[6] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[11] * (xtau[21] - x[6]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[11] * (0.35e1 * x[6] / (0.30e1 * x[6] + 0.5e0 * x[7] + 0.5e0) - x[6]);
  xdot[7] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[11] * (xtau[22] - x[7]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[11] * (x[7] / (0.30e1 * x[6] + 0.5e0 * x[7] + 0.5e0) - x[7]);
  xdot[8] = 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[11] * (xtau[23] - x[8]) - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[11], -0.1e1 / 0.3e1) * (x[8] - alpha[10]) / 0.100e3;
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
plhs[0] = mxCreateDoubleMatrix(1,(mwSize)9,mxREAL);                                                                                       
/* get a pointer to the real data in the output matrix */                                                                                 
double *xdot = mxGetPr(plhs[0]);                                                                                                          
DDErightHandSide(x, xtau, alpha, p, xdot);                                                                                                   
}                                                                                                                                         

