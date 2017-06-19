%% mex-Wrapper for DDE

fileID = fopen('DDE_wrapped.c','w');
fprintf(fileID,'#include "mex.h"\n');
addID = fopen('DDE.c','r');
text = fread(addID,'*char');
fclose(addID);
fprintf(fileID,text);
clear text addID
fprintf(fileID,'\n\n');
fprintf(fileID,'/* The gateway function */\n');
fprintf(fileID,'void mexFunction (int nlhs, mxArray *plhs[],\n');
fprintf(fileID,'int nrhs, const mxArray *prhs[])\n');
fprintf(fileID,'{\n');
fprintf(fileID,'if(nrhs!=4) {\n');
fprintf(fileID,'mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nrhs","4 inputs required (state vector, matrix with delayed states and vector with parameters).");\n');
fprintf(fileID,'}\n');
fprintf(fileID,'if(nlhs==0) {\n');
fprintf(fileID,'mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nlhs","Please define an output!");\n');
fprintf(fileID,'}\n');
fprintf(fileID,'if(nlhs!=1) {\n');
fprintf(fileID,'mexErrMsgIdAndTxt("MyToolbox:DDEWrapped:nlhs","One output required.");\n');
fprintf(fileID,'}\n');
fprintf(fileID,'double *x = mxGetPr(prhs[0]);\n');
fprintf(fileID,'double *xtau = mxGetPr(prhs[1]);\n');
fprintf(fileID,'double *alpha = mxGetPr(prhs[2]);\n');
fprintf(fileID,'double *p = mxGetPr(prhs[3]);\n');
fprintf(fileID,'/* create the output matrix */\n');
fprintf(fileID,'plhs[0] = mxCreateDoubleMatrix(1,(mwSize)%d,mxREAL);\n',xnum);       %%% %d vorher 12
fprintf(fileID,'/* get a pointer to the real data in the output matrix */\n');
fprintf(fileID,'double *xdot = mxGetPr(plhs[0]);\n');
fprintf(fileID,'DDErightHandSide(x, xtau, alpha, xdot);\n');
fprintf(fileID,'}\n');