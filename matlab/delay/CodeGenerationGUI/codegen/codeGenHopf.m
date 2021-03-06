%> @file codeHopf.m
%> @brief if hopf manifold was chosen, this code will be added to the maple
%> code

%% Code generation Hopf

fileID = fopen(cg_name,'a');
fprintf(fileID,'# define normal vector system\n');
fprintf(fileID,'AugSys:=AugSys2:-SdDelayBif:-HopfNV:-CreateHopfNVSys(Sys,[ '); %%%%
for i=1:anum-1
fprintf(fileID,'alpha%d, ',i);
end
fprintf(fileID,'alpha%d]):-getSys();\n',anum);
fprintf(fileID,'\n');
fprintf(fileID,'# pick relevant equations of Normal Vectors System\n');
fprintf(fileID,'manifoldEq:=[]:\n');
fprintf(fileID,'for i from 1 by 1 to %d do\n',(3*xnum+2));
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
fprintf(fileID,'Procedure4CodeGen:=proc(x, alpha,p, omega, w1, w2 )\n');
fprintf(fileID,'m;\n');
fprintf(fileID,'end proc;\n');
fprintf(fileID,'ManifoldEquation:=subs([m=manifoldEq],eval(Procedure4CodeGen));\n');
fprintf(fileID,'# generate C code\n');
fprintf(fileID,'CodeGeneration:-C(ManifoldEquation,returnvariablename="residuum",defaulttype=numeric,output="HopfManifold.c",deducetypes=false);\n');
fprintf(fileID,'\n');
% NV
fprintf(fileID,'# pick relevant equations of Normal Vectors System\n');
fprintf(fileID,'manifoldEq:=[]:\n');
fprintf(fileID,'for i from %d by 1 to nops(AugSys["Equations"]) do\n',(3*xnum+3));
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
fprintf(fileID,'Procedure4CodeGen:=proc(x, alpha, omega,p, w1, w2, v1, v2, g1, g2, u, r)\n');
fprintf(fileID,'m;\n');
fprintf(fileID,'end proc;\n');
fprintf(fileID,'ManifoldEquation:=subs([m=manifoldEq],eval(Procedure4CodeGen));\n');
fprintf(fileID,'# generate C code\n');
fprintf(fileID,'CodeGeneration:-C(ManifoldEquation,returnvariablename="residuum",defaulttype=numeric,output="HopfNV.c",deducetypes=false);\n');
fprintf(fileID,'\n');
fprintf(fileID,'\n\n\n');
