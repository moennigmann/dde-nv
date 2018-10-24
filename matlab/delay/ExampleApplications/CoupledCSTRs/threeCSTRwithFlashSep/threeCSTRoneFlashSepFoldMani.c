#include "mex.h"
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double *w1,
  double residuum[25])
{
  residuum[0] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (0.1e1 - x[0]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 * x[9] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[0]) - 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0];
  residuum[1] = -0.5e1 / 0.18e2 * alpha[3] / alpha[5] * x[1] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[10] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[1]) + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * x[0] - 0.2500e4 * exp(-0.7216742844e4 / x[2]) * x[1];
  residuum[2] = 0.5e1 / 0.18e2 * alpha[3] / alpha[5] * (alpha[2] - x[2]) + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (x[11] - x[2]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * x[0] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) * (x[2] - alpha[4]) / 0.10e2;
  residuum[3] = 0.5e1 / 0.18e2 * alpha[7] / alpha[9] * (0.1e1 - x[3]) + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * (x[0] - x[3]) - 0.2770e4 * exp(-0.6013952370e4 / x[5]) * x[3];
  residuum[4] = -0.5e1 / 0.18e2 * alpha[7] / alpha[9] * x[4] + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * (x[1] - x[4]) + 0.2770e4 * exp(-0.6013952370e4 / x[5]) * x[3] - 0.2500e4 * exp(-0.7216742844e4 / x[5]) * x[4];
  residuum[5] = 0.5e1 / 0.18e2 * alpha[7] / alpha[9] * (alpha[6] - x[5]) + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * (x[2] - x[5]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[5]) * x[3] + 0.4166666675e8 * exp(-0.7216742844e4 / x[5]) * x[4] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[9], -0.1e1 / 0.3e1) * (x[5] - alpha[8]) / 0.10e2;
  residuum[6] = 0.5e1 / 0.18e2 * alpha[11] / alpha[13] * (0.1e1 - x[6]) + 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] * (x[3] - x[6]) - 0.2770e4 * exp(-0.6013952370e4 / x[8]) * x[6];
  residuum[7] = -0.5e1 / 0.18e2 * alpha[11] / alpha[13] * x[7] + 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] * (x[4] - x[7]) + 0.2770e4 * exp(-0.6013952370e4 / x[8]) * x[6] - 0.2500e4 * exp(-0.7216742844e4 / x[8]) * x[7];
  residuum[8] = 0.5e1 / 0.18e2 * alpha[11] / alpha[13] * (alpha[10] - x[8]) + 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] * (x[5] - x[8]) + 0.3957142861e8 * exp(-0.6013952370e4 / x[8]) * x[6] + 0.4166666675e8 * exp(-0.7216742844e4 / x[8]) * x[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[13], -0.1e1 / 0.3e1) * (x[8] - alpha[12]) / 0.10e2;
  residuum[9] = 0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] * (x[6] - x[9]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[15] * (0.35e1 * x[9] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[9]);
  residuum[10] = 0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] * (x[7] - x[10]) - 0.5e1 / 0.18e2 * alpha[0] / alpha[15] * (x[10] / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - x[10]);
  residuum[11] = 0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] * (x[8] - x[11]) - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[15], -0.1e1 / 0.3e1) * (x[11] - alpha[14]) / 0.10e2;
  residuum[12] = (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2770e4 * exp(-0.6013952370e4 / x[2])) * w1[0] + 0.2770e4 * exp(-0.6013952370e4 / x[2]) * w1[1] + 0.3957142861e8 * exp(-0.6013952370e4 / x[2]) * w1[2] + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * w1[3];
  residuum[13] = (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] - 0.2500e4 * exp(-0.7216742844e4 / x[2])) * w1[1] + 0.4166666675e8 * exp(-0.7216742844e4 / x[2]) * w1[2] + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * w1[4];
  residuum[14] = -0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] * w1[0] + (0.1665864806e8 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] - 0.1804185711e8 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1]) * w1[1] + (-0.5e1 / 0.18e2 * alpha[3] / alpha[5] - 0.5e1 / 0.18e2 * alpha[0] / alpha[5] + 0.2379806869e12 * pow(x[2], -0.2e1) * exp(-0.6013952370e4 / x[2]) * x[0] + 0.3006976191e12 * pow(x[2], -0.2e1) * exp(-0.7216742844e4 / x[2]) * x[1] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[5], -0.1e1 / 0.3e1) / 0.10e2) * w1[2] + 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] * w1[5];
  residuum[15] = (-0.5e1 / 0.18e2 * alpha[7] / alpha[9] - 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] - 0.2770e4 * exp(-0.6013952370e4 / x[5])) * w1[3] + 0.2770e4 * exp(-0.6013952370e4 / x[5]) * w1[4] + 0.3957142861e8 * exp(-0.6013952370e4 / x[5]) * w1[5] + 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] * w1[6];
  residuum[16] = (-0.5e1 / 0.18e2 * alpha[7] / alpha[9] - 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] - 0.2500e4 * exp(-0.7216742844e4 / x[5])) * w1[4] + 0.4166666675e8 * exp(-0.7216742844e4 / x[5]) * w1[5] + 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] * w1[7];
  residuum[17] = -0.1665864806e8 * pow(x[5], -0.2e1) * exp(-0.6013952370e4 / x[5]) * x[3] * w1[3] + (0.1665864806e8 * pow(x[5], -0.2e1) * exp(-0.6013952370e4 / x[5]) * x[3] - 0.1804185711e8 * pow(x[5], -0.2e1) * exp(-0.7216742844e4 / x[5]) * x[4]) * w1[4] + (-0.5e1 / 0.18e2 * alpha[7] / alpha[9] - 0.1000e4 * (alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[9] + 0.2379806869e12 * pow(x[5], -0.2e1) * exp(-0.6013952370e4 / x[5]) * x[3] + 0.3006976191e12 * pow(x[5], -0.2e1) * exp(-0.7216742844e4 / x[5]) * x[4] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[9], -0.1e1 / 0.3e1) / 0.10e2) * w1[5] + 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] * w1[8];
  residuum[18] = (-0.5e1 / 0.18e2 * alpha[11] / alpha[13] - 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] - 0.2770e4 * exp(-0.6013952370e4 / x[8])) * w1[6] + 0.2770e4 * exp(-0.6013952370e4 / x[8]) * w1[7] + 0.3957142861e8 * exp(-0.6013952370e4 / x[8]) * w1[8] + 0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] * w1[9];
  residuum[19] = (-0.5e1 / 0.18e2 * alpha[11] / alpha[13] - 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] - 0.2500e4 * exp(-0.7216742844e4 / x[8])) * w1[7] + 0.4166666675e8 * exp(-0.7216742844e4 / x[8]) * w1[8] + 0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] * w1[10];
  residuum[20] = -0.1665864806e8 * pow(x[8], -0.2e1) * exp(-0.6013952370e4 / x[8]) * x[6] * w1[6] + (0.1665864806e8 * pow(x[8], -0.2e1) * exp(-0.6013952370e4 / x[8]) * x[6] - 0.1804185711e8 * pow(x[8], -0.2e1) * exp(-0.7216742844e4 / x[8]) * x[7]) * w1[7] + (-0.5e1 / 0.18e2 * alpha[11] / alpha[13] - 0.1000e4 * (alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[13] + 0.2379806869e12 * pow(x[8], -0.2e1) * exp(-0.6013952370e4 / x[8]) * x[6] + 0.3006976191e12 * pow(x[8], -0.2e1) * exp(-0.7216742844e4 / x[8]) * x[7] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[13], -0.1e1 / 0.3e1) / 0.10e2) * w1[8] + 0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] * w1[11];
  residuum[21] = (-0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] - 0.5e1 / 0.18e2 * alpha[0] / alpha[15] * (0.35e1 / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - 0.1050e2 * x[9] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1) - 0.1e1)) * w1[9] + 0.8333333333e0 * alpha[0] / alpha[15] * x[10] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1) * w1[10] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.35e1 / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - 0.1050e2 * x[9] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1)) * w1[0] - 0.8333333333e0 * alpha[0] / alpha[5] * x[10] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1) * w1[1];
  residuum[22] = 0.4861111111e0 * alpha[0] / alpha[15] * x[9] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1) * w1[9] + (-0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] - 0.5e1 / 0.18e2 * alpha[0] / alpha[15] * (0.1e1 / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - 0.5e0 * x[10] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1) - 0.1e1)) * w1[10] - 0.4861111111e0 * alpha[0] / alpha[5] * x[9] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1) * w1[0] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * (0.1e1 / (0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0) - 0.5e0 * x[10] * pow(0.30e1 * x[9] + 0.5e0 * x[10] + 0.5e0, -0.2e1)) * w1[1];
  residuum[23] = (-0.1000e4 * (alpha[11] / 0.3600e4 + alpha[7] / 0.3600e4 + alpha[3] / 0.3600e4 + alpha[0] / 0.3600e4) / alpha[15] - alpha[1] * pow(0.1000e4, 0.1e1 / 0.3e1) * pow(alpha[15], -0.1e1 / 0.3e1) / 0.10e2) * w1[11] + 0.5e1 / 0.18e2 * alpha[0] / alpha[5] * w1[2];
  residuum[24] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) + pow(w1[6], 0.2e1) + pow(w1[7], 0.2e1) + pow(w1[8], 0.2e1) + pow(w1[9], 0.2e1) + pow(w1[10], 0.2e1) + pow(w1[11], 0.2e1) - 0.1e1;
}

/* The gateway function */                                                                                               
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])                                             
{                                                                                                                        
                                                                                                                         
 /* check for proper number of arguments */                                                                              
 if(nrhs!=4) {                                                                                                           
     mexErrMsgIdAndTxt("MyToolbox:threeCSTRoneFlashSepHopfManiPop:nrhs","4 inputs required (some of them are vectors).");
 }                                                                                                                       
 if(nlhs==0) {                                                                                                           
     mexErrMsgIdAndTxt("MyToolbox:threeCSTRoneFlashSepHopfManiPop:nlhs","Please define an output!");                     
 }                                                                                                                       
 if(nlhs!=1) {                                                                                                           
     mexErrMsgIdAndTxt("MyToolbox:threeCSTRoneFlashSepHopfManiPop:nlhs","One output required.");                         
 }                                                                                                                       
                                                                                                                         
 /* get the values of the inputs                                                                                         
  * scalars are fetched using mxGetScalar and assigned to a type double,                                                 
  * vectors are fetched using mxGetPr and assigned to type double pointer                                                
  */                                                                                                                     
                                                                                                                         
 double *xPointer = mxGetPr(prhs[0]);                                                                                    
 double *alphaPointer = mxGetPr(prhs[1]);                                                                                      
 double *pPointer = mxGetPr(prhs[2]);                                                                                     
 double *w1Pointer = mxGetPr(prhs[3]);                                                                                     
                                                                                                                         
                                                                                                                         
 /* create the output matrix */                                                                                          
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)                                                                               
25                                                                                                                       
,mxREAL);                                                                                                                
                                                                                                                         
 /* get a pointer to the real data in the output matrix */                                                               
 double *residuumPointer = mxGetPr(plhs[0]);                                                                             
                                                                                                                         
 /* call the computational routine */                                                                                    
 ManifoldEquation(xPointer, alphaPointer, w1Pointer, residuumPointer);                                   
}                                                                                                                        

