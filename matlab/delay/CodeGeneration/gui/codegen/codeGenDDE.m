%> @file codeGenDDE.m
%> @brief adds the dde-code to the maple code

%% Code generation - get dde equations

fileID = fopen(cg_name,'a');
fprintf(fileID,'#get dde equations\n');
fprintf(fileID,'dderhs:=[]; # create space for dde right hand side\n');
fprintf(fileID,'for ii from 1 by 1 to nops(Sys["DynVars"]) do\n');
fprintf(fileID,'ithRhs:=rhs(Sys["ODEs"][ii]): # get one entry and replace...\n');
fprintf(fileID,'for jj from 1 by 1 to nops(Sys["DynVars"]) do\n');
fprintf(fileID,'ithRhs:=subs(Sys["DynVars"][jj]=x[jj],ithRhs): # ...states, ...\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'for jj from 1 by 1 to nops(Sys["DelVars"]) do\n');
fprintf(fileID,'ithRhs:=subs(Sys["DelVars"][jj]=xtau[jj],ithRhs); # ...delayed states,...\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'for jj from 1 by 1 to nops(Sys["Parameters"]) do\n');
fprintf(fileID,'ithRhs:=subs(lhs(Sys["Parameters"][jj])=alpha[jj],ithRhs): # and parameters...\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'dderhs:=[op(dderhs),ithRhs]; # and concatenate it with the existing entries\n');
fprintf(fileID,'end do:\n');
fprintf(fileID,'# create frame for code generation\n');
fprintf(fileID,'Procedure4CodeGen:=proc(x,xtau,alpha)\n');
fprintf(fileID,'m;\n');
fprintf(fileID,'end proc;\n');
fprintf(fileID,'DDErightHandSide:=subs([m=dderhs],eval(Procedure4CodeGen));\n');
fprintf(fileID,'# generate C code\n');
fprintf(fileID,'CodeGeneration:-C(DDErightHandSide,returnvariablename="xdot",defaulttype=numeric,output="DDE.c",deducetypes=false);');