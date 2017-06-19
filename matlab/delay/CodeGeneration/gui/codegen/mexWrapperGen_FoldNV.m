%% Mex Wrapper - Fold NV
% generates c-code to be compilable in mex-file for matlab
% fold


fileID = fopen('fold_NV.c','w');
fprintf(fileID,'#include "mex.h"\n');
addID = fopen('FoldNV.c','r');
text = fread(addID,'*char');
fclose(addID);
fprintf(fileID,text);
clear text addID
fprintf(fileID,'\n\n');
fprintf(fileID,'/* The gateway function */\n');
fprintf(fileID,'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) \n');
fprintf(fileID,'{\n\n');
fprintf(fileID,' /* check for proper number of arguments */\n');
fprintf(fileID,' if(nrhs!=5) {  \n'); % hier 5?
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldNV:nrhs","5 inputs required (some of them are vectors).");\n'); % liste
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
fprintf(fileID,' double *rPointer = mxGetPr(prhs[4]);\n\n\n');
fprintf(fileID,' /* create the output matrix */\n');
fprintf(fileID,' plhs[0] = mxCreateDoubleMatrix(1,(mwSize)%d,mxREAL);\n\n',(2*xnum+anum+1)); %%% vorher 7
fprintf(fileID,' /* get a pointer to the real data in the output matrix */\n');
fprintf(fileID,' double *residuumPointer = mxGetPr(plhs[0]);\n\n');
fprintf(fileID,' /* call the computational routine */\n');
fprintf(fileID,' ManifoldEquation(xPointer,alphaPointer, pPointer, wPointer,vPointer,g1,uPointer,rPointer, residuumPointer);\n');
fprintf(fileID,'}');