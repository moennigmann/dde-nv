%% Code generation initialisation

fileID = fopen('codegen.txt','w');
fprintf(fileID,'#initialisation\n');
fprintf(fileID,'restart; \n');
fprintf(fileID,'#define paths\n');
fprintf(fileID,'_ModulesDirectory:="');
fprintf(fileID,'%s";\n',way);
fprintf(fileID,'#load modules \n cat(_ModulesDirectory, "/Aux/Aux.mpl");\n read(cat(_ModulesDirectory, "/Aux/Aux.mpl")): \n read(cat(_ModulesDirectory, "/AugSys2/AugSys2.mpl")): \n'); 
fprintf(fileID,'\n');
fprintf(fileID,'#define system\n');