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
  double residuum[15])
{
  residuum[0] = -0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[0] + 0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[1] - 0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[3] + 0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[4] - 0.1000e4 * x[2] * v[0] - (-0.1000e4 * p[0] + 0.5000e4 * x[2]) * v[1] + 0.2000e1 * (0.2e1 * x[2] + 0.1e1) * x[0] * v[2] - 0.10e1 * v[0] + 0.2e1 * g1 * w[0];
  residuum[1] = -0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[0] - 0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[1] - 0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[3] - 0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[4] - (0.1000e4 * p[0] - 0.5000e4 * x[2]) * v[0] - 0.1000e4 * x[2] * v[1] + 0.2000e1 * (0.2e1 * x[2] + 0.1e1) * x[1] * v[2] - 0.10e1 * v[1] + 0.2e1 * g1 * w[1];
  residuum[2] = -(0.1000e4 * x[0] - 0.5000e4 * x[1]) * v[0] - (0.5000e4 * x[0] + 0.1000e4 * x[1]) * v[1] - (-0.1000e1 - 0.2000e1 * pow(x[0], 0.2e1) - 0.2000e1 * pow(x[1], 0.2e1)) * v[2] - 0.10e1 * v[2] + 0.2e1 * g1 * w[2];
  residuum[3] = -0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[0] + 0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[1] - 0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[3] + 0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[4] - 0.1000e4 * x[5] * v[3] - (-0.1000e4 * p[0] + 0.5000e4 * x[5]) * v[4] + 0.2000e1 * (0.2e1 * x[5] + 0.1e1) * x[3] * v[5] - 0.10e1 * v[3] + 0.2e1 * g1 * w[3];
  residuum[4] = -0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[0] - 0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[1] - 0.1359140914e1 * sin(0.1000e4 * p[0] + 0.1e1) * v[3] - 0.1359140914e1 * cos(0.1000e4 * p[0] + 0.1e1) * v[4] - (0.1000e4 * p[0] - 0.5000e4 * x[5]) * v[3] - 0.1000e4 * x[5] * v[4] + 0.2000e1 * (0.2e1 * x[5] + 0.1e1) * x[4] * v[5] - 0.10e1 * v[4] + 0.2e1 * g1 * w[4];
  residuum[5] = -(0.1000e4 * x[3] - 0.5000e4 * x[4]) * v[3] - (0.5000e4 * x[3] + 0.1000e4 * x[4]) * v[4] - (-0.1000e1 - 0.2000e1 * pow(x[3], 0.2e1) - 0.2000e1 * pow(x[4], 0.2e1)) * v[5] - 0.10e1 * v[5] + 0.2e1 * g1 * w[5];
  residuum[6] = 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[0] - 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[1] + 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[3] - 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[4] + 0.1000e4 * x[2] * u[0] + (-0.1000e4 * p[0] + 0.5000e4 * x[2]) * u[1] - 0.2000e1 * (0.2e1 * x[2] + 0.1e1) * x[0] * u[2] - 0.1000e4 * v[0] * w[2] - 0.5000e4 * v[1] * w[2] - v[2] * ((-0.4000e1 * x[2] - 0.2000e1) * w[0] - 0.4000e1 * x[0] * w[2]);
  residuum[7] = 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[0] + 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[1] + 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[3] + 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[4] + (0.1000e4 * p[0] - 0.5000e4 * x[2]) * u[0] + 0.1000e4 * x[2] * u[1] - 0.2000e1 * (0.2e1 * x[2] + 0.1e1) * x[1] * u[2] + 0.5000e4 * v[0] * w[2] - 0.1000e4 * v[1] * w[2] - v[2] * ((-0.4000e1 * x[2] - 0.2000e1) * w[1] - 0.4000e1 * x[1] * w[2]);
  residuum[8] = (0.1000e4 * x[0] - 0.5000e4 * x[1]) * u[0] + (0.5000e4 * x[0] + 0.1000e4 * x[1]) * u[1] + (-0.1000e1 - 0.2000e1 * pow(x[0], 0.2e1) - 0.2000e1 * pow(x[1], 0.2e1)) * u[2] - v[0] * (0.1000e4 * w[0] - 0.5000e4 * w[1]) - v[1] * (0.5000e4 * w[0] + 0.1000e4 * w[1]) - v[2] * (-0.4000e1 * x[0] * w[0] - 0.4000e1 * x[1] * w[1]);
  residuum[9] = 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[0] - 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[1] + 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[3] - 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[4] + 0.1000e4 * x[5] * u[3] + (-0.1000e4 * p[0] + 0.5000e4 * x[5]) * u[4] - 0.2000e1 * (0.2e1 * x[5] + 0.1e1) * x[3] * u[5] - 0.1000e4 * v[3] * w[5] - 0.5000e4 * v[4] * w[5] - v[5] * ((-0.4000e1 * x[5] - 0.2000e1) * w[3] - 0.4000e1 * x[3] * w[5]);
  residuum[10] = 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[0] + 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[1] + 0.5000e0 * sin(0.1000e4 * p[0] + 0.1e1) * u[3] + 0.5000e0 * cos(0.1000e4 * p[0] + 0.1e1) * u[4] + (0.1000e4 * p[0] - 0.5000e4 * x[5]) * u[3] + 0.1000e4 * x[5] * u[4] - 0.2000e1 * (0.2e1 * x[5] + 0.1e1) * x[4] * u[5] + 0.5000e4 * v[3] * w[5] - 0.1000e4 * v[4] * w[5] - v[5] * ((-0.4000e1 * x[5] - 0.2000e1) * w[4] - 0.4000e1 * x[4] * w[5]);
  residuum[11] = (0.1000e4 * x[3] - 0.5000e4 * x[4]) * u[3] + (0.5000e4 * x[3] + 0.1000e4 * x[4]) * u[4] + (-0.1000e1 - 0.2000e1 * pow(x[3], 0.2e1) - 0.2000e1 * pow(x[4], 0.2e1)) * u[5] - v[3] * (0.1000e4 * w[3] - 0.5000e4 * w[4]) - v[4] * (0.5000e4 * w[3] + 0.1000e4 * w[4]) - v[5] * (-0.4000e1 * x[3] * w[3] - 0.4000e1 * x[4] * w[4]);
  residuum[12] = 0.1000e1 * u[2] - r[0];
  residuum[13] = 0.1000e1 * u[5] - r[1];
  residuum[14] = pow(r[0], 0.2e1) + pow(r[1], 0.2e1) - 0.1e1;
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
 plhs[0] = mxCreateDoubleMatrix(1,(mwSize)15,mxREAL);                                                                 
                                                                                                                     
 /* get a pointer to the real data in the output matrix */                                                           
 double *residuumPointer = mxGetPr(plhs[0]);                                                                         
                                                                                                                     
 /* call the computational routine */                                                                                
 ManifoldEquation(xPointer,alphaPointer, pPointer, wPointer,vPointer,g1,uPointer,rPointer, residuumPointer);                   
}
