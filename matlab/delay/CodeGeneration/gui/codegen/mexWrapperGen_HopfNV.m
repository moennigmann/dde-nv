%% c-code fuer matrlab holen -> ziel: sc3hopfNV.c
% generates c-code to be compilable in mex-file for matlab


fileID = fopen('hopf_NV.c','w');
fprintf(fileID,'#include "mex.h"\n');
addID = fopen('HopfNV.c','r');
text = fread(addID,'*char');
fclose(addID);
fprintf(fileID,text);
clear text addID
fprintf(fileID,'\n\n');
fprintf(fileID,'/* The gateway function */\n');
fprintf(fileID,'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\n');
fprintf(fileID,'{\n\n');
fprintf(fileID,'  /* check for proper number of arguments */\n');
fprintf(fileID,' if(nrhs!=12) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nrhs","12 inputs required (some of them are vectors).");\n');
fprintf(fileID,' }\n');
fprintf(fileID,' if(nlhs==0) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nlhs","Please define an output!");\n');
fprintf(fileID,' }\n');
fprintf(fileID,' if(nlhs!=1) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:GenEigManiPop:nlhs","One output required.");\n');
fprintf(fileID,' }\n\n');
fprintf(fileID,' /* get the values of the inputs\n');
fprintf(fileID,'  * scalars are fetched using mxGetScalar and assigned to a type double,\n');
fprintf(fileID,'  * vectors are fetched using mxGetPr and assigned to type double pointer\n');
fprintf(fileID,'  */\n\n');
fprintf(fileID,' double *xPointer = mxGetPr(prhs[0]);\n');
fprintf(fileID,' double *alphaPointer = mxGetPr(prhs[1]);\n');
fprintf(fileID,' double *pPointer = mxGetPr(prhs[2]);\n');
fprintf(fileID,' double omega = mxGetScalar(prhs[3]);\n');
fprintf(fileID,' double *w1Pointer = mxGetPr(prhs[4]);\n');
fprintf(fileID,' double *w2Pointer = mxGetPr(prhs[5]);\n');
fprintf(fileID,' double *v1Pointer = mxGetPr(prhs[6]);\n');
fprintf(fileID,' double *v2Pointer = mxGetPr(prhs[7]);\n');
fprintf(fileID,' double g1 = mxGetScalar(prhs[8]);\n');
fprintf(fileID,' double g2 = mxGetScalar(prhs[9]);\n');
fprintf(fileID,' double *uPointer = mxGetPr(prhs[10]);\n');
fprintf(fileID,' double *rPointer = mxGetPr(prhs[11]);\n\n\n');
fprintf(fileID,' /* create the output matrix */\n');
fprintf(fileID,' plhs[0] = mxCreateDoubleMatrix(1,(mwSize)%d,mxREAL);\n\n',(6*xnum+anum+4)); % vorher 47
fprintf(fileID,' /* get a pointer to the real data in the output matrix */\n');
fprintf(fileID,' double *residuumPointer = mxGetPr(plhs[0]);\n\n');
fprintf(fileID,' /* call the computational routine */\n');
fprintf(fileID,' ManifoldEquation(xPointer,alphaPointer,omega,w1Pointer,w2Pointer,v1Pointer,\n');
fprintf(fileID,'         v2Pointer,g1,g2,uPointer,rPointer,residuumPointer);\n');
fprintf(fileID,'}\n');