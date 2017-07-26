%> @file codeGenInit.m
%> @brief initialisation code is generated, including inserting the path to
%> the maple files

%% Code generation initialisation

fileID = fopen(cg_name,'w');
fprintf(fileID,'#initialisation\n');
fprintf(fileID,'restart; \n');
fprintf(fileID,'#define paths\n');
fprintf(fileID,'_ModulesDirectory:="');
fprintf(fileID,'%s";\n',way);
fprintf(fileID,'#load modules \n cat(_ModulesDirectory, "/Aux/Aux.mpl");\n read(cat(_ModulesDirectory, "/Aux/Aux.mpl")): \n read(cat(_ModulesDirectory, "/AugSys2/AugSys2.mpl")): \n'); 
fprintf(fileID,'\n');
fprintf(fileID,'#define system\n');