%% Mex Wrapper - modFold NV
% generates c-code to be compilable in mex-file for matlab
% -> mod fold?


fileID = fopen('modfold_NV.c','w');
fprintf(fileID,'#include "mex.h"\n');
addID = fopen('NV.c','r');
text = fread(addID,'*char');
fclose(addID);
fprintf(fileID,text);
clear text addID
fprintf(fileID,'\n\n');
fprintf(fileID,'/* The gateway function */\n');
fprintf(fileID,'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) \n');
fprintf(fileID,'{\n\n');
fprintf(fileID,' /* check for proper number of arguments */\n');
fprintf(fileID,' if(nrhs!=8) {  \n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldNV:nrhs","8 inputs required (some of them are vectors).");\n');
fprintf(fileID,' }\n');
fprintf(fileID,' if(nlhs==0) { \n');
fprintf(fileID,'    mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldNV:nlhs","Please define an output!");\n');
fprintf(fileID,' }\n');
fprintf(fileID,' if(nlhs!=1) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldNV:nlhs","One output required.");\n');
fprintf(fileID,' }\n\n');
fprintf(fileID,' /* get the values of the inputs \n');
fprintf(fileID,'  * scalars are fetched using mxGetScalar and assigned to a type double,\n');
fprintf(fileID,'  * vectors are fetched using mxGetPr and assigned to type double pointer\n');
fprintf(fileID,'  */\n\n');
fprintf(fileID,' double *xPointer = mxGetPr(prhs[0]);\n');
fprintf(fileID,' double *alphaPointer = mxGetPr(prhs[1]);\n');
fprintf(fileID,' double *pPointer = mxGetPr(prhs[2]);\n');
fprintf(fileID,' double *wPointer = mxGetPr(prhs[3]);\n');
fprintf(fileID,' double *vPointer = mxGetPr(prhs[4]);\n');
fprintf(fileID,' double g1 = mxGetScalar(prhs[5]);\n');
fprintf(fileID,' double *uPointer = mxGetPr(prhs[6]);\n');
fprintf(fileID,' double *rPointer = mxGetPr(prhs[7]);\n\n\n');
fprintf(fileID,' /* create the output matrix */\n');
fprintf(fileID,' plhs[0] = mxCreateDoubleMatrix(1,(mwSize)%d,mxREAL);\n\n',(4*xnum+2*anum+1)); % war 7
fprintf(fileID,' /* get a pointer to the real data in the output matrix */\n');
fprintf(fileID,' double *residuumPointer = mxGetPr(plhs[0]);\n\n');
fprintf(fileID,' /* call the computational routine */\n');
fprintf(fileID,' ManifoldEquation(xPointer,alphaPointer, pPointer, wPointer,vPointer,g1,uPointer,rPointer, residuumPointer);\n');
fprintf(fileID,'}');