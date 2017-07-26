%> @file mexWrapperGen_Modfold.m
%> @brief if modfold manifold was chosen, this code will be added to the compiled 
%> maple code, to make it compilable to matlabcode

%% Mex-Wrapper - modFold
% generates c-code to be compilable in mex-file for matlab

fileID = fopen('modfold_Manifold.c','w');
fprintf(fileID,'#include "mex.h"\n');
addID = fopen('Manifold.c','r');
text = fread(addID,'*char');
fclose(addID);
fprintf(fileID,text);
clear text addID
fprintf(fileID,'\n\n');
fprintf(fileID,'/* The gateway function */\n');
fprintf(fileID,'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\n');
fprintf(fileID,'{\n\n');
fprintf(fileID,' /* check for proper number of arguments */\n');
fprintf(fileID,' if(nrhs!=4) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldMani:nrhs","4 inputs required (some of them are vectors).");\n');
fprintf(fileID,' }\n');
fprintf(fileID,' if(nlhs==0) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldMani:nlhs","Please define an output!");\n');
fprintf(fileID,' }\n');
fprintf(fileID,' if(nlhs!=1) {\n');
fprintf(fileID,'     mexErrMsgIdAndTxt("MyToolbox:DDENLPmodfoldMani:nlhs","One output required.");\n');
fprintf(fileID,' }\n\n');
fprintf(fileID,' /* get the values of the inputs \n');
fprintf(fileID,'  * scalars are fetched using mxGetScalar and assigned to a type double,\n');
fprintf(fileID,'  * vectors are fetched using mxGetPr and assigned to type double pointer\n');
fprintf(fileID,'  */\n\n');
fprintf(fileID,' double *xPointer = mxGetPr(prhs[0]);\n');
fprintf(fileID,' double *alphaPointer = mxGetPr(prhs[1]);\n');
fprintf(fileID,' double *pPointer = mxGetPr(prhs[2]);\n');
fprintf(fileID,' double *wPointer = mxGetPr(prhs[3]);\n\n\n');
fprintf(fileID,' /* create the output matrix */\n');
fprintf(fileID,' plhs[0] = mxCreateDoubleMatrix(1,(mwSize)\n');
fprintf(fileID,'%d\n',(2*xnum+1));      
fprintf(fileID,',mxREAL);\n\n');
fprintf(fileID,' /* get a pointer to the real data in the output matrix */\n');
fprintf(fileID,' double *residuumPointer = mxGetPr(plhs[0]);\n\n');
fprintf(fileID,' /* call the computational routine */\n');
fprintf(fileID,' ManifoldEquation(xPointer, alphaPointer, wPointer, residuumPointer);\n'); 
fprintf(fileID,'}');