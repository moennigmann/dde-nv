#include "mex.h"
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double *w1,
  double *r,
  double residuum[16])
{
  residuum[ 0] = 0.5e1 / 0.18e2 / alpha[5] * (0.35e1 * x[9] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[0]) * w1[0] + 0.5e1 / 0.18e2 / alpha[5] * (x[10] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[1]) * w1[1] + 0.5e1 / 0.18e2 / alpha[5] * (x[11] - x[2]) * w1[2] + 0.5e1 / 0.18e2 / alpha[9] * (x[0] - x[3]) * w1[3] + 0.5e1 / 0.18e2 / alpha[9] * (x[1] - x[4]) * w1[4] + 0.5e1 / 0.18e2 / alpha[9] * (x[2] - x[5]) * w1[5] + 0.5e1 / 0.18e2 / alpha[13] * (x[3] - x[6]) * w1[6] + 0.5e1 / 0.18e2 / alpha[13] * (x[4] - x[7]) * w1[7] + 0.5e1 / 0.18e2 / alpha[13] * (x[5] - x[8]) * w1[8] + (0.5e1 / 0.18e2 / alpha[15] * (x[6] - x[9]) - 0.5e1 / 0.18e2 / alpha[15] * (0.35e1 * x[9] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[9])) * w1[9] + (0.5e1 / 0.18e2 / alpha[15] * (x[7] - x[10]) - 0.5e1 / 0.18e2 / alpha[15] * (x[10] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[10])) * w1[10] + 0.5e1 / 0.18e2 / alpha[15] * (x[8] - x[11]) * w1[11] - r[0];
  residuum[ 1] = -pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) * w1[2] / 0.10e2 - pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[9], -0.1e1 / 0.3e1) * (x[5] - alpha[8]) * w1[5] / 0.10e2 - pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[13], -0.1e1 / 0.3e1) * (x[8] - alpha[12]) * w1[8] / 0.10e2 - pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[15], -0.1e1 / 0.3e1) * (x[11] - alpha[14]) * w1[11] / 0.10e2 - r[1];
  residuum[ 2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * w1[2] - r[2];
  residuum[ 3] = 0.5e1 / 0.18e2 / alpha[5] * (0.1e1 - x[0]) * w1[0] - 0.5e1 / 0.18e2 / alpha[5] * x[1] * w1[1] + 0.5e1 / 0.18e2 / alpha[5] * (alpha[2] - x[2]) * w1[2] + 0.5e1 / 0.18e2 / alpha[9] * (x[0] - x[3]) * w1[3] + 0.5e1 / 0.18e2 / alpha[9] * (x[1] - x[4]) * w1[4] + 0.5e1 / 0.18e2 / alpha[9] * (x[2] - x[5]) * w1[5] + 0.5e1 / 0.18e2 / alpha[13] * (x[3] - x[6]) * w1[6] + 0.5e1 / 0.18e2 / alpha[13] * (x[4] - x[7]) * w1[7] + 0.5e1 / 0.18e2 / alpha[13] * (x[5] - x[8]) * w1[8] + 0.5e1 / 0.18e2 / alpha[15] * (x[6] - x[9]) * w1[9] + 0.5e1 / 0.18e2 / alpha[15] * (x[7] - x[10]) * w1[10] + 0.5e1 / 0.18e2 / alpha[15] * (x[8] - x[11]) * w1[11] - r[3];
  residuum[ 4] = alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * w1[2] / 0.10e2 - r[4];
  residuum[ 5] = (-0.5e1 / 0.18e2 * alpha[3] * pow(alpha[5], -0.2e1) * (0.1e1 - x[0]) - 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[5], -0.2e1) * (0.35e1 * x[9] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[0])) * w1[0] + (0.5e1 / 0.18e2 * alpha[3] * pow(alpha[5], -0.2e1) * x[1] - 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[5], -0.2e1) * (x[10] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[1])) * w1[1] + (-0.5e1 / 0.18e2 * alpha[3] * pow(alpha[5], -0.2e1) * (alpha[2] - x[2]) - 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[5], -0.2e1) * (x[11] - x[2]) + alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.4e1 / 0.3e1) * (x[2] - alpha[4]) / 0.30e2) * w1[2] - r[5];
  residuum[ 6] = 0.5e1 / 0.18e2 * alpha[7] / alpha[9] * w1[5] - r[6];
  residuum[ 7] = 0.5e1 / 0.18e2 / alpha[9] * (0.1e1 - x[3]) * w1[3] - 0.5e1 / 0.18e2 / alpha[9] * x[4] * w1[4] + 0.5e1 / 0.18e2 / alpha[9] * (alpha[6] - x[5]) * w1[5] + 0.5e1 / 0.18e2 / alpha[13] * (x[3] - x[6]) * w1[6] + 0.5e1 / 0.18e2 / alpha[13] * (x[4] - x[7]) * w1[7] + 0.5e1 / 0.18e2 / alpha[13] * (x[5] - x[8]) * w1[8] + 0.5e1 / 0.18e2 / alpha[15] * (x[6] - x[9]) * w1[9] + 0.5e1 / 0.18e2 / alpha[15] * (x[7] - x[10]) * w1[10] + 0.5e1 / 0.18e2 / alpha[15] * (x[8] - x[11]) * w1[11] - r[7];
  residuum[ 8] = alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[9], -0.1e1 / 0.3e1) * w1[5] / 0.10e2 - r[8];
  residuum[ 9] = (-0.5e1 / 0.18e2 * alpha[7] * pow(alpha[9], -0.2e1) * (0.1e1 - x[3]) - 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[9], -0.2e1) * (x[0] - x[3])) * w1[3] + (0.5e1 / 0.18e2 * alpha[7] * pow(alpha[9], -0.2e1) * x[4] - 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[9], -0.2e1) * (x[1] - x[4])) * w1[4] + (-0.5e1 / 0.18e2 * alpha[7] * pow(alpha[9], -0.2e1) * (alpha[6] - x[5]) - 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[9], -0.2e1) * (x[2] - x[5]) + alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[9], -0.4e1 / 0.3e1) * (x[5] - alpha[8]) / 0.30e2) * w1[5] - r[9];
  residuum[10] = 0.5e1 / 0.18e2 * alpha[11] / alpha[13] * w1[8] - r[10];
  residuum[11] = 0.5e1 / 0.18e2 / alpha[13] * (0.1e1 - x[6]) * w1[6] - 0.5e1 / 0.18e2 / alpha[13] * x[7] * w1[7] + 0.5e1 / 0.18e2 / alpha[13] * (alpha[10] - x[8]) * w1[8] + 0.5e1 / 0.18e2 / alpha[15] * (x[6] - x[9]) * w1[9] + 0.5e1 / 0.18e2 / alpha[15] * (x[7] - x[10]) * w1[10] + 0.5e1 / 0.18e2 / alpha[15] * (x[8] - x[11]) * w1[11] - r[11];
  residuum[12] = alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[13], -0.1e1 / 0.3e1) * w1[8] / 0.10e2 - r[12];
  residuum[13] = (-0.5e1 / 0.18e2 * alpha[11] * pow(alpha[13], -0.2e1) * (0.1e1 - x[6]) - 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[13], -0.2e1) * (x[3] - x[6])) * w1[6] + (0.5e1 / 0.18e2 * alpha[11] * pow(alpha[13], -0.2e1) * x[7] - 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[13], -0.2e1) * (x[4] - x[7])) * w1[7] + (-0.5e1 / 0.18e2 * alpha[11] * pow(alpha[13], -0.2e1) * (alpha[10] - x[8]) - 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[13], -0.2e1) * (x[5] - x[8]) + alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[13], -0.4e1 / 0.3e1) * (x[8] - alpha[12]) / 0.30e2) * w1[8] - r[13];
  residuum[14] = alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[15], -0.1e1 / 0.3e1) * w1[11] / 0.10e2 - r[14];
  residuum[15] = (-0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[15], -0.2e1) * (x[6] - x[9]) + 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[15], -0.2e1) * (0.35e1 * x[9] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[9])) * w1[9] + (-0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[15], -0.2e1) * (x[7] - x[10]) + 0.5e1 / 0.18e2 * alpha[0] * pow(alpha[15], -0.2e1) * (x[10] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[10])) * w1[10] + (-0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) * pow(alpha[15], -0.2e1) * (x[8] - x[11]) + alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[15], -0.4e1 / 0.3e1) * (x[11] - alpha[14]) / 0.30e2) * w1[11] - r[15];
}

/* The gateway function */                                                                                                       
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])                                                     
{                                                                                                                                
                                                                                                                                 
 /* check for proper number of arguments */                                                                                      
 if(nrhs!=5) {                                                                                                                  
     mexErrMsgIdAndTxt("MyToolbox:threeCSTRoneFlashSepHopfNVPop:nrhs","5 inputs required (some of them are vectors).");         
 }                                                                                                                               
 if(nlhs==0) {                                                                                                                   
    mexErrMsgIdAndTxt("MyToolbox:threeCSTRoneFlashSepHopfNVPop:nlhs","Please define an output!");                                
 }                                                                                                                               
 if(nlhs!=1) {                                                                                                                   
     mexErrMsgIdAndTxt("MyToolbox:threeCSTRoneFlashSepHopfNVPop:nlhs","One output required.");                                   
 }                                                                                                                               
                                                                                                                                 
 /* get the values of the inputs                                                                                                 
  * scalars are fetched using mxGetScalar and assigned to a type double,                                                         
  * vectors are fetched using mxGetPr and assigned to type double pointer                                                        
  */                                                                                                                             
                                                                                                                                 
 double *xPointer = mxGetPr(prhs[0]);                                                                                            
 double *alphaPointer = mxGetPr(prhs[1]);                                                                                       
 double *pPointer = mxGetPr(prhs[2]);                                                                                             
 double *w1Pointer = mxGetPr(prhs[3]);                                                                                           
 double *rPointer = mxGetPr(prhs[4]);                                                                                           
                                                                                                                                 
                                                                                                                                 
 /* create the output matrix */                                                                                                  
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)16,mxREAL);                                                                            
                                                                                                                                 
 /* get a pointer to the real data in the output matrix */                                                                       
 double *residuumPointer = mxGetPr(plhs[0]);                                                                                     
                                                                                                                                 
 /* call the computational routine */                                                                                            
 ManifoldEquation(xPointer,alphaPointer,w1Pointer,rPointer,residuumPointer);
}                                                                                                                                

