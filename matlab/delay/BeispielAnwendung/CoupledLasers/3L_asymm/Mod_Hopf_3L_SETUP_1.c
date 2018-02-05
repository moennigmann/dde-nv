#include "mex.h"
#include <math.h>

void ManifoldEquation (
  double *x,
  double *alpha,
  double *p,
  double omega,
  double *w1,
  double *w2,
  double residuum[29])
{
  residuum[0] = 0.1000e4 * p[0] * x[1] + 0.5000e3 * x[2] * x[0] - 0.20000e4 * x[2] * x[1] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * x[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * x[1];
  residuum[1] = -0.1000e4 * p[0] * x[0] + 0.20000e4 * x[2] * x[0] + 0.5000e3 * x[2] * x[1] - 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * x[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * x[1];
  residuum[2] = 0.5000e1 * alpha[0] - 0.5000e1 * x[2] - 0.5000e1 * (x[2] + 0.1e1) * (pow(x[0], 0.2e1) + pow(x[1], 0.2e1));
  residuum[3] = 0.1000e4 * p[0] * x[4] + 0.5000e3 * x[5] * x[3] - 0.20000e4 * x[5] * x[4] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * x[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * x[1];
  residuum[4] = -0.1000e4 * p[0] * x[3] + 0.20000e4 * x[5] * x[3] + 0.5000e3 * x[5] * x[4] - 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * x[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * x[1];
  residuum[5] = 0.5000e1 * alpha[1] - 0.5000e1 * x[5] - 0.5000e1 * (x[5] + 0.1e1) * (pow(x[3], 0.2e1) + pow(x[4], 0.2e1));
  residuum[6] = 0.1000e4 * p[0] * x[7] + 0.5000e3 * x[8] * x[6] - 0.20000e4 * x[8] * x[7] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * x[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * x[1];
  residuum[7] = -0.1000e4 * p[0] * x[6] + 0.20000e4 * x[8] * x[6] + 0.5000e3 * x[8] * x[7] - 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * x[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * x[1];
  residuum[8] = 0.5000e1 * alpha[2] - 0.5000e1 * x[8] - 0.5000e1 * (x[8] + 0.1e1) * (pow(x[6], 0.2e1) + pow(x[7], 0.2e1));
  residuum[9] = -0.1105170918e1 * cos(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.1105170918e1 * sin(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[1]) - 0.5000e3 * x[2] * w1[0] - (0.1000e4 * p[0] - 0.20000e4 * x[2]) * w1[1] - (0.5000e3 * x[0] - 0.20000e4 * x[1]) * w1[2] - omega * w2[0] - 0.10e1 * w1[0];
  residuum[10] = -0.1105170918e1 * cos(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.1105170918e1 * sin(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[1]) - (-0.1000e4 * p[0] + 0.20000e4 * x[2]) * w1[0] - 0.5000e3 * x[2] * w1[1] - (0.20000e4 * x[0] + 0.5000e3 * x[1]) * w1[2] - omega * w2[1] - 0.10e1 * w1[1];
  residuum[11] = 0.10000e2 * (x[2] + 0.1e1) * x[0] * w1[0] + 0.10000e2 * (x[2] + 0.1e1) * x[1] * w1[1] - (-0.5000e1 - 0.5000e1 * pow(x[0], 0.2e1) - 0.5000e1 * pow(x[1], 0.2e1)) * w1[2] - omega * w2[2] - 0.10e1 * w1[2];
  residuum[12] = -0.1105170918e1 * cos(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.1105170918e1 * sin(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[1]) - 0.5000e3 * x[5] * w1[3] - (0.1000e4 * p[0] - 0.20000e4 * x[5]) * w1[4] - (0.5000e3 * x[3] - 0.20000e4 * x[4]) * w1[5] - omega * w2[3] - 0.10e1 * w1[3];
  residuum[13] = -0.1105170918e1 * cos(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.1105170918e1 * sin(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[1]) - (-0.1000e4 * p[0] + 0.20000e4 * x[5]) * w1[3] - 0.5000e3 * x[5] * w1[4] - (0.20000e4 * x[3] + 0.5000e3 * x[4]) * w1[5] - omega * w2[4] - 0.10e1 * w1[4];
  residuum[14] = 0.10000e2 * (x[5] + 0.1e1) * x[3] * w1[3] + 0.10000e2 * (x[5] + 0.1e1) * x[4] * w1[4] - (-0.5000e1 - 0.5000e1 * pow(x[3], 0.2e1) - 0.5000e1 * pow(x[4], 0.2e1)) * w1[5] - omega * w2[5] - 0.10e1 * w1[5];
  residuum[15] = -0.1105170918e1 * cos(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.1105170918e1 * sin(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[1]) - 0.5000e3 * x[8] * w1[6] - (0.1000e4 * p[0] - 0.20000e4 * x[8]) * w1[7] - (0.5000e3 * x[6] - 0.20000e4 * x[7]) * w1[8] - omega * w2[6] - 0.10e1 * w1[6];
  residuum[16] = -0.1105170918e1 * cos(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.1105170918e1 * sin(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[1]) - (-0.1000e4 * p[0] + 0.20000e4 * x[8]) * w1[6] - 0.5000e3 * x[8] * w1[7] - (0.20000e4 * x[6] + 0.5000e3 * x[7]) * w1[8] - omega * w2[7] - 0.10e1 * w1[7];
  residuum[17] = 0.10000e2 * (x[8] + 0.1e1) * x[6] * w1[6] + 0.10000e2 * (x[8] + 0.1e1) * x[7] * w1[7] - (-0.5000e1 - 0.5000e1 * pow(x[6], 0.2e1) - 0.5000e1 * pow(x[7], 0.2e1)) * w1[8] - omega * w2[8] - 0.10e1 * w1[8];
  residuum[18] = -0.1105170918e1 * cos(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[1]) + 0.1105170918e1 * sin(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.5000e3 * x[2] * w2[0] - (0.1000e4 * p[0] - 0.20000e4 * x[2]) * w2[1] - (0.5000e3 * x[0] - 0.20000e4 * x[1]) * w2[2] + omega * w1[0] - 0.10e1 * w2[0];
  residuum[19] = -0.1105170918e1 * cos(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[1]) + 0.1105170918e1 * sin(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[1]) - (-0.1000e4 * p[0] + 0.20000e4 * x[2]) * w2[0] - 0.5000e3 * x[2] * w2[1] - (0.20000e4 * x[0] + 0.5000e3 * x[1]) * w2[2] + omega * w1[1] - 0.10e1 * w2[1];
  residuum[20] = 0.10000e2 * (x[2] + 0.1e1) * x[0] * w2[0] + 0.10000e2 * (x[2] + 0.1e1) * x[1] * w2[1] - (-0.5000e1 - 0.5000e1 * pow(x[0], 0.2e1) - 0.5000e1 * pow(x[1], 0.2e1)) * w2[2] + omega * w1[2] - 0.10e1 * w2[2];
  residuum[21] = -0.1105170918e1 * cos(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[1]) + 0.1105170918e1 * sin(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.5000e3 * x[5] * w2[3] - (0.1000e4 * p[0] - 0.20000e4 * x[5]) * w2[4] - (0.5000e3 * x[3] - 0.20000e4 * x[4]) * w2[5] + omega * w1[3] - 0.10e1 * w2[3];
  residuum[22] = -0.1105170918e1 * cos(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[1]) + 0.1105170918e1 * sin(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[1]) - (-0.1000e4 * p[0] + 0.20000e4 * x[5]) * w2[3] - 0.5000e3 * x[5] * w2[4] - (0.20000e4 * x[3] + 0.5000e3 * x[4]) * w2[5] + omega * w1[4] - 0.10e1 * w2[4];
  residuum[23] = 0.10000e2 * (x[5] + 0.1e1) * x[3] * w2[3] + 0.10000e2 * (x[5] + 0.1e1) * x[4] * w2[4] - (-0.5000e1 - 0.5000e1 * pow(x[3], 0.2e1) - 0.5000e1 * pow(x[4], 0.2e1)) * w2[5] + omega * w1[5] - 0.10e1 * w2[5];
  residuum[24] = -0.1105170918e1 * cos(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[1]) + 0.1105170918e1 * sin(omega / 0.10e2) * (0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[1]) - 0.5000e3 * x[8] * w2[6] - (0.1000e4 * p[0] - 0.20000e4 * x[8]) * w2[7] - (0.5000e3 * x[6] - 0.20000e4 * x[7]) * w2[8] + omega * w1[6] - 0.10e1 * w2[6];
  residuum[25] = -0.1105170918e1 * cos(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w2[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w2[1]) + 0.1105170918e1 * sin(omega / 0.10e2) * (-0.5000e1 * sin(0.2e1 + 0.100e3 * p[0]) * w1[0] + 0.5000e1 * cos(0.2e1 + 0.100e3 * p[0]) * w1[1]) - (-0.1000e4 * p[0] + 0.20000e4 * x[8]) * w2[6] - 0.5000e3 * x[8] * w2[7] - (0.20000e4 * x[6] + 0.5000e3 * x[7]) * w2[8] + omega * w1[7] - 0.10e1 * w2[7];
  residuum[26] = 0.10000e2 * (x[8] + 0.1e1) * x[6] * w2[6] + 0.10000e2 * (x[8] + 0.1e1) * x[7] * w2[7] - (-0.5000e1 - 0.5000e1 * pow(x[6], 0.2e1) - 0.5000e1 * pow(x[7], 0.2e1)) * w2[8] + omega * w1[8] - 0.10e1 * w2[8];
  residuum[27] = pow(w1[0], 0.2e1) + pow(w1[1], 0.2e1) + pow(w1[2], 0.2e1) + pow(w1[3], 0.2e1) + pow(w1[4], 0.2e1) + pow(w1[5], 0.2e1) + pow(w1[6], 0.2e1) + pow(w1[7], 0.2e1) + pow(w1[8], 0.2e1) + pow(w2[0], 0.2e1) + pow(w2[1], 0.2e1) + pow(w2[2], 0.2e1) + pow(w2[3], 0.2e1) + pow(w2[4], 0.2e1) + pow(w2[5], 0.2e1) + pow(w2[6], 0.2e1) + pow(w2[7], 0.2e1) + pow(w2[8], 0.2e1) - 0.1e1;
  residuum[28] = w1[0] * w2[0] + w1[1] * w2[1] + w1[2] * w2[2] + w1[3] * w2[3] + w1[4] * w2[4] + w1[5] * w2[5] + w1[6] * w2[6] + w1[7] * w2[7] + w1[8] * w2[8];
}


/* The gateway function */                                                                                             
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])                                           
{                                                                                                                      
                                                                                                                       
 /* check for proper number of arguments */                                                                            
 if(nrhs!=6) {                                                                                                         
     mexErrMsgIdAndTxt("MyToolbox:oneCSTRoneFlashSepHopfManiPop:nrhs","6 inputs required (some of them are vectors).");
 }                                                                                                                     
 if(nlhs==0) {                                                                                                         
     mexErrMsgIdAndTxt("MyToolbox:oneCSTRoneFlashSepHopfManiPop:nlhs","Please define an output!");                     
 }                                                                                                                     
 if(nlhs!=1) {                                                                                                         
     mexErrMsgIdAndTxt("MyToolbox:oneCSTRoneFlashSepHopfManiPop:nlhs","One output required.");                         
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
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)                                                                             
29                                                                                                                     
,mxREAL);                                                                                                              
                                                                                                                       
 /* get a pointer to the real data in the output matrix */                                                             
 double *residuumPointer = mxGetPr(plhs[0]);                                                                           
                                                                                                                       
 /* call the computational routine */                                                                                  
 ManifoldEquation(xPointer, alphaPointer, pPointer, omega,w1Pointer,w2Pointer, residuumPointer);                                 
}  