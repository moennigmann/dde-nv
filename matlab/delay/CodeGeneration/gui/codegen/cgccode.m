%% Code generation c-code Fold NV

fileID = fopen('codegen.txt','a');
fprintf(fileID,'#pick relevant equations of Normal Vextor System\n');
fprintf(fileID,'manifoldEq:=[]:\n');
fprintf(fileID,'for i from %d by 1 to nops(AugSys["Equations"]) do\n',(2*xnum+2));
fprintf(fileID,'ithRhs:=rhs(AugSys["Equations"][i]);\n');
fprintf(fileID,'for jj from 1 by 1 to nops(Sys["DynVars"]) do\n');
fprintf(fileID,'ithRhs:=subs(Sys["DynVars"][jj]=x[jj],ithRhs): # ...states, ...\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'for jj from 1 by 1 to nops(Sys["Parameters"]) do\n');
fprintf(fileID,'ithRhs:=subs(lhs(Sys["Parameters"][jj])=alpha[jj],ithRhs): # and parameters...\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'manifoldEq:=[op(manifoldEq),ithRhs]:\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'# create frame for code generation\n');
fprintf(fileID,'Procedure4CodeGen:=proc(x, alpha,  w1, r)\n');
fprintf(fileID,'m;\n');
fprintf(fileID,'end proc;\n');
fprintf(fileID,'ManifoldEquation:=subs([m=manifoldEq],eval(Procedure4CodeGen));\n');
fprintf(fileID,'# generate C code\n');
fprintf(fileID,'CodeGeneration:-C(ManifoldEquation,returnvariablename="residuum",defaulttype=numeric,output="FoldNV.c",deducetypes=false);\n');
fprintf(fileID,'\n\n\n');