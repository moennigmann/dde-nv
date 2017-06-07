%% Code generation define system
% comment line 50-52 and 59-66 to make code general

fileID = fopen('codegen.txt','a');
fprintf(fileID,'Sys["DynVars"]:=[');
%% print states for DynVars
% xnames(i,1) contains counted states
if xnum>0           % check for input
    for i=1:xnum-1
        fprintf(fileID,'%s,',cell2mat(xnames(i,1)));
    end
    fprintf(fileID,'%s]; \n',cell2mat(xnames(xnum,1)));
else 
    fprintf(fileID,'];\n');
end
%
fprintf(fileID,'Sys["Parameters"]:=[');
%% print parameters (counted) with their value
% alphavec(i,1) contains counted parameter alphavec(i,3) contains the value
if anum>0           % check for input
    for i=1:anum-1
        fprintf(fileID,'%s=%d,',cell2mat(alphavec(i,1)),cell2mat(alphavec(i,3)));
    end
    fprintf(fileID,'%s=%d]; \n',cell2mat(alphavec(anum,1)),cell2mat(alphavec(anum,3)));
else
    fprintf(fileID,'];\n');
end
%
fprintf(fileID,'Sys["AEs"]:=[];\n');
fprintf(fileID,'Sys["ODEs"]:=[');
%% print system equations
% xnames(i,1) contains conted states xdot(i,2) contains equations
if xnum>0           % check for input
    for i=1:xnum-1
        fprintf(fileID,'`%s''` = %s,\n',cell2mat(xnames(i,1)),cell2mat(xdot(i,2)));
    end
    fprintf(fileID,'`%s''` = %s]; \n',cell2mat(xnames(xnum,1)),cell2mat(xdot(xnum,2)));
else 
    fprintf(fileID,'];\n');
end
%
fprintf(fileID,'Sys["DelVars"]:=[');
%% print delay variables for states
% xnames(i,1) contains counted states tauname(j) contains counted delays
if xnum>0           % check for input
    for j=1:delnum-1
        for i=1:xnum
        fprintf(fileID,'%s%s, ',cell2mat(xnames(i,1)),cell2mat(tauname(j)));
        end
    end
    for i=1:xnum-1
        fprintf(fileID,'%s%s, ',cell2mat(xnames(i,1)),cell2mat(tauname(delnum)));
    end
fprintf(fileID,'%s%s]; \n',cell2mat(xnames(xnum,1)),cell2mat(tauname(delnum)));
%
else 
    fprintf(fileID,'];\n');
end
fprintf(fileID,'Sys["AlgVars"]:=[];\n');
%% print explicit AEs
% xnames(i,2) contains names of the states 
% xnames(i,1) contains counted states
% alphavec(i,2) contains names of the parameter
% alphavec(i,1) contains counted parameter
% aes(i) contains additional equations (eg: parameter in the ODEs)
% tauname(j) contains counted delays
fprintf(fileID,'Sys["ExplicitAEs"]:=[');
% for i=1:14
%     fprintf(fileID,'%s = %s, ',cell2mat(p(i,1)),cell2mat(p(i,2)));
% end
% print xname = xi
for i=1:xnum
   fprintf(fileID,'%s=%s, ',cell2mat(xnames(i,2)),cell2mat(xnames(i,1)));
end
% print alphaname = alphai
for i=1:anum
    fprintf(fileID,'%s=%s, ',cell2mat(alphavec(i,2)),cell2mat(alphavec(i,1)));
end
% fprintf(fileID,'alphaA=p1, alphaB=p2, alphaC=p3, ');
% fprintf(fileID,'EA=p4, EB=p5, ');
% fprintf(fileID,'Cp=p6, dHA=p7, dHB=p8, BA=p9, BB=p10, ');
% fprintf(fileID,'R=p11, kA=p12, kB=p13, Fp=p14, ');
% fprintf(fileID,'cA1tau1=x1tau1, cB1tau1=x2tau1, T1tau1=x3tau1, ');
% fprintf(fileID,'cAFtau1=x4tau1, cBFtau1=x5tau1, TFtau1=x6tau1, ');
% fprintf(fileID,'cA1tau2=x1tau2, cB1tau2=x2tau2, T1tau2=x3tau2, ');
% fprintf(fileID,'cAFtau2=x4tau2, cBFtau2=x5tau2, TFtau2=x6tau2];\n');
% print additional explicit algebraic equations
if l            % check whether used or not
for i = 1:l
    fprintf(fileID,'%s, ',cell2mat(aes(i)));
end
end
% print xnametauj = xitauj
for j=1:delnum-1
    for i=1:xnum
        fprintf(fileID,'%s%s=%s%s, ',cell2mat(xnames(i,2)),cell2mat(tauname(j)),cell2mat(xnames(i,1)),cell2mat(tauname(j)));
    end
end
    for i=1:xnum-1
        fprintf(fileID,'%s%s=%s%s, ',cell2mat(xnames(i,2)),cell2mat(tauname(delnum)),cell2mat(xnames(i,1)),cell2mat(tauname(delnum)));
    end
fprintf(fileID,'%s%s=%s%s];\n',cell2mat(xnames(xnum,2)),cell2mat(tauname(delnum)),cell2mat(xnames(xnum,1)),cell2mat(tauname(delnum)));
% ];\n');
%
%% print expressions for delays
fprintf(fileID,'Sys["Delays"]:=[');
for i=1:(delnum-1)
    fprintf(fileID,'tau[%d]=%s, ',i,cell2mat(del(i)));
end
fprintf(fileID,'tau[%d]=%s];\n ',delnum,cell2mat(del(delnum)));
fprintf(fileID,'\n');
%
%% print insert fixed parameters
fprintf(fileID,'#insert fixed parameters\n#\n');
fprintf(fileID,'Sys:=Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(Sys);\n');
fprintf(fileID,'# look for errors\n');
fprintf(fileID,'Aux:-SystemClasses:-listOfErrorsInDDESys(Sys, strict);\n');
% fprintf(fileID,'# define normal vector system\n');
% fprintf(fileID,'AugSys:=AugSys2:-SdDelayBif:-FoldNV:-CreateFoldNVSys(Sys,[ '); %%%%
% for i=1:anum-1
% fprintf(fileID,'alpha%d, ',i);
% end
% fprintf(fileID,'alpha%d]):-getSys();\n',anum);
% fprintf(fileID,'\n');
%