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
  double residuum[22])
{
  residuum[0] = -0.5525854590e1 * cos(0.100e3 * p[0] + 0.2e1) * v[0] + 0.5525854590e1 * sin(0.100e3 * p[0] + 0.2e1) * v[1] - 0.5525854590e1 * cos(0.100e3 * p[0] + 0.2e1) * v[3] + 0.5525854590e1 * sin(0.100e3 * p[0] + 0.2e1) * v[4] - 0.5525854590e1 * cos(0.100e3 * p[0] + 0.2e1) * v[6] + 0.5525854590e1 * sin(0.100e3 * p[0] + 0.2e1) * v[7] - 0.1000e4 * x[2] * v[0] - (-0.1000e4 * p[0] + 0.4000e4 * x[2]) * v[1] + 0.10000e2 * (0.2e1 * x[2] + 0.1e1) * x[0] * v[2] - 0.1e1 * v[0] + 0.2e1 * g1 * w[0];
  residuum[1] = -0.5525854590e1 * sin(0.100e3 * p[0] + 0.2e1) * v[0] - 0.5525854590e1 * cos(0.100e3 * p[0] + 0.2e1) * v[1] - 0.5525854590e1 * sin(0.100e3 * p[0] + 0.2e1) * v[3] - 0.5525854590e1 * cos(0.100e3 * p[0] + 0.2e1) * v[4] - 0.5525854590e1 * sin(0.100e3 * p[0] + 0.2e1) * v[6] - 0.5525854590e1 * cos(0.100e3 * p[0] + 0.2e1) * v[7] - (0.1000e4 * p[0] - 0.4000e4 * x[2]) * v[0] - 0.1000e4 * x[2] * v[1] + 0.10000e2 * (0.2e1 * x[2] + 0.1e1) * x[1] * v[2] - 0.1e1 * v[1] + 0.2e1 * g1 * w[1];
  residuum[2] = -(0.1000e4 * x[0] - 0.4000e4 * x[1]) * v[0] - (0.4000e4 * x[0] + 0.1000e4 * x[1]) * v[1] - (-0.5000e1 - 0.10000e2 * pow(x[0], 0.2e1) - 0.10000e2 * pow(x[1], 0.2e1)) * v[2] - 0.1e1 * v[2] + 0.2e1 * g1 * w[2];
  residuum[3] = -0.1000e4 * x[5] * v[3] - (-0.1000e4 * p[0] + 0.4000e4 * x[5]) * v[4] + 0.10000e2 * (0.2e1 * x[5] + 0.1e1) * x[3] * v[5] - 0.1e1 * v[3] + 0.2e1 * g1 * w[3];
  residuum[4] = -(0.1000e4 * p[0] - 0.4000e4 * x[5]) * v[3] - 0.1000e4 * x[5] * v[4] + 0.10000e2 * (0.2e1 * x[5] + 0.1e1) * x[4] * v[5] - 0.1e1 * v[4] + 0.2e1 * g1 * w[4];
  residuum[5] = -(0.1000e4 * x[3] - 0.4000e4 * x[4]) * v[3] - (0.4000e4 * x[3] + 0.1000e4 * x[4]) * v[4] - (-0.5000e1 - 0.10000e2 * pow(x[3], 0.2e1) - 0.10000e2 * pow(x[4], 0.2e1)) * v[5] - 0.1e1 * v[5] + 0.2e1 * g1 * w[5];
  residuum[6] = -0.1000e4 * x[8] * v[6] - (-0.1000e4 * p[0] + 0.4000e4 * x[8]) * v[7] + 0.10000e2 * (0.2e1 * x[8] + 0.1e1) * x[6] * v[8] - 0.1e1 * v[6] + 0.2e1 * g1 * w[6];
  residuum[7] = -(0.1000e4 * p[0] - 0.4000e4 * x[8]) * v[6] - 0.1000e4 * x[8] * v[7] + 0.10000e2 * (0.2e1 * x[8] + 0.1e1) * x[7] * v[8] - 0.1e1 * v[7] + 0.2e1 * g1 * w[7];
  residuum[8] = -(0.1000e4 * x[6] - 0.4000e4 * x[7]) * v[6] - (0.4000e4 * x[6] + 0.1000e4 * x[7]) * v[7] - (-0.5000e1 - 0.10000e2 * pow(x[6], 0.2e1) - 0.10000e2 * pow(x[7], 0.2e1)) * v[8] - 0.1e1 * v[8] + 0.2e1 * g1 * w[8];
  residuum[9] = 0.5000e1 * cos(0.100e3 * p[0] + 0.2e1) * u[0] - 0.5000e1 * sin(0.100e3 * p[0] + 0.2e1) * u[1] + 0.5000e1 * cos(0.100e3 * p[0] + 0.2e1) * u[3] - 0.5000e1 * sin(0.100e3 * p[0] + 0.2e1) * u[4] + 0.5000e1 * cos(0.100e3 * p[0] + 0.2e1) * u[6] - 0.5000e1 * sin(0.100e3 * p[0] + 0.2e1) * u[7] + 0.1000e4 * x[2] * u[0] + (-0.1000e4 * p[0] + 0.4000e4 * x[2]) * u[1] - 0.10000e2 * (0.2e1 * x[2] + 0.1e1) * x[0] * u[2] - 0.1000e4 * v[0] * w[2] - 0.4000e4 * v[1] * w[2] - v[2] * ((-0.20000e2 * x[2] - 0.10000e2) * w[0] - 0.20000e2 * x[0] * w[2]);
  residuum[10] = 0.5000e1 * sin(0.100e3 * p[0] + 0.2e1) * u[0] + 0.5000e1 * cos(0.100e3 * p[0] + 0.2e1) * u[1] + 0.5000e1 * sin(0.100e3 * p[0] + 0.2e1) * u[3] + 0.5000e1 * cos(0.100e3 * p[0] + 0.2e1) * u[4] + 0.5000e1 * sin(0.100e3 * p[0] + 0.2e1) * u[6] + 0.5000e1 * cos(0.100e3 * p[0] + 0.2e1) * u[7] + (0.1000e4 * p[0] - 0.4000e4 * x[2]) * u[0] + 0.1000e4 * x[2] * u[1] - 0.10000e2 * (0.2e1 * x[2] + 0.1e1) * x[1] * u[2] + 0.4000e4 * v[0] * w[2] - 0.1000e4 * v[1] * w[2] - v[2] * ((-0.20000e2 * x[2] - 0.10000e2) * w[1] - 0.20000e2 * x[1] * w[2]);
  residuum[11] = (0.1000e4 * x[0] - 0.4000e4 * x[1]) * u[0] + (0.4000e4 * x[0] + 0.1000e4 * x[1]) * u[1] + (-0.5000e1 - 0.10000e2 * pow(x[0], 0.2e1) - 0.10000e2 * pow(x[1], 0.2e1)) * u[2] - v[0] * (0.1000e4 * w[0] - 0.4000e4 * w[1]) - v[1] * (0.4000e4 * w[0] + 0.1000e4 * w[1]) - v[2] * (-0.20000e2 * x[0] * w[0] - 0.20000e2 * x[1] * w[1]);
  residuum[12] = 0.1000e4 * x[5] * u[3] + (-0.1000e4 * p[0] + 0.4000e4 * x[5]) * u[4] - 0.10000e2 * (0.2e1 * x[5] + 0.1e1) * x[3] * u[5] - 0.1000e4 * v[3] * w[5] - 0.4000e4 * v[4] * w[5] - v[5] * ((-0.20000e2 * x[5] - 0.10000e2) * w[3] - 0.20000e2 * x[3] * w[5]);
  residuum[13] = (0.1000e4 * p[0] - 0.4000e4 * x[5]) * u[3] + 0.1000e4 * x[5] * u[4] - 0.10000e2 * (0.2e1 * x[5] + 0.1e1) * x[4] * u[5] + 0.4000e4 * v[3] * w[5] - 0.1000e4 * v[4] * w[5] - v[5] * ((-0.20000e2 * x[5] - 0.10000e2) * w[4] - 0.20000e2 * x[4] * w[5]);
  residuum[14] = (0.1000e4 * x[3] - 0.4000e4 * x[4]) * u[3] + (0.4000e4 * x[3] + 0.1000e4 * x[4]) * u[4] + (-0.5000e1 - 0.10000e2 * pow(x[3], 0.2e1) - 0.10000e2 * pow(x[4], 0.2e1)) * u[5] - v[3] * (0.1000e4 * w[3] - 0.4000e4 * w[4]) - v[4] * (0.4000e4 * w[3] + 0.1000e4 * w[4]) - v[5] * (-0.20000e2 * x[3] * w[3] - 0.20000e2 * x[4] * w[4]);
  residuum[15] = 0.1000e4 * x[8] * u[6] + (-0.1000e4 * p[0] + 0.4000e4 * x[8]) * u[7] - 0.10000e2 * (0.2e1 * x[8] + 0.1e1) * x[6] * u[8] - 0.1000e4 * v[6] * w[8] - 0.4000e4 * v[7] * w[8] - v[8] * ((-0.20000e2 * x[8] - 0.10000e2) * w[6] - 0.20000e2 * x[6] * w[8]);
  residuum[16] = (0.1000e4 * p[0] - 0.4000e4 * x[8]) * u[6] + 0.1000e4 * x[8] * u[7] - 0.10000e2 * (0.2e1 * x[8] + 0.1e1) * x[7] * u[8] + 0.4000e4 * v[6] * w[8] - 0.1000e4 * v[7] * w[8] - v[8] * ((-0.20000e2 * x[8] - 0.10000e2) * w[7] - 0.20000e2 * x[7] * w[8]);
  residuum[17] = (0.1000e4 * x[6] - 0.4000e4 * x[7]) * u[6] + (0.4000e4 * x[6] + 0.1000e4 * x[7]) * u[7] + (-0.5000e1 - 0.10000e2 * pow(x[6], 0.2e1) - 0.10000e2 * pow(x[7], 0.2e1)) * u[8] - v[6] * (0.1000e4 * w[6] - 0.4000e4 * w[7]) - v[7] * (0.4000e4 * w[6] + 0.1000e4 * w[7]) - v[8] * (-0.20000e2 * x[6] * w[6] - 0.20000e2 * x[7] * w[7]);
  residuum[18] = 0.5000e1 * u[2] - r[0];
  residuum[19] = 0.5000e1 * u[5] - r[1];
  residuum[20] = 0.5000e1 * u[8] - r[2];
  residuum[21] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) + pow(r[2], 0.2e1) - 0.1e1;
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
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)22,mxREAL);                                                                 
                                                                                                                     
 /* get a pointer to the real data in the output matrix */                                                           
 double *residuumPointer = mxGetPr(plhs[0]);                                                                         
                                                                                                                     
 /* call the computational routine */                                                                                
 ManifoldEquation(xPointer,alphaPointer, pPointer, wPointer,vPointer,g1,uPointer,rPointer, residuumPointer);                   
}
