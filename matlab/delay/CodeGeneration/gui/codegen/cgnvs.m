%% Code generation insert fixed parameters and define normal vector system

fileID = fopen('codegen.txt','a');
fprintf(fileID,'#insert fixed parameters\n#\n');
fprintf(fileID,'Sys:=Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(Sys);\n');
fprintf(fileID,'# look for errors\n');
fprintf(fileID,'Aux:-SystemClasses:-listOfErrorsInDDESys(Sys, strict);\n');
fprintf(fileID,'# define normal vector system\n');
fprintf(fileID,'AugSys:=AugSys2:-SdDelayBif:-FoldNV:-CreateFoldNVSys(Sys,[ ');
for i=1:anum-1
fprintf(fileID,'alpha%d, ',i);
end
fprintf(fileID,'alpha%d]):-getSys();\n',anum);
fprintf(fileID,'\n');