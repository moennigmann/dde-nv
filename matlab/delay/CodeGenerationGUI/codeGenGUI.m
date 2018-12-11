%> @file codeGenGUI.m
%> @brief Programm to control input dialogs (gui), code generation and
%> compiling

%% Program to open input dialogs (gui), code generation, code manipulation and compiling
% add folders gui and codegen to matlab

%% start gui

systemInput;
% systemInput_limit;

%% generate code for maple

codeGeneration;
%%
% clear aes alphavec ans answer del fileID i j l tauname way x xdot xnames maxrealpart

system_command = ['maple ' cg_name];

system(system_command);

clear p_name system_order
%% code manipulation to make it compilable for matlab

if choiceMani(1) % Fold
    mexWrapperGen_Fold;
    mexWrapperGen_FoldNV;
end
if choiceMani(2) % modFold
    mexWrapperGen_Modfold;
    mexWrapperGen_ModfoldNV;
end
if choiceMani(3) % Hopf
    mexWrapperGen_HopfMani;
    mexWrapperGen_HopfNV;
end
if choiceMani(4) % modHopf
    mexWrapperGen_ModHopf;
    mexWrapperGen_ModHopfNV;
end
mexWrapperGen_DDE;
%% Compile for matlab

% replace some code
cfiles = dir('*.c');
for cfile = cfiles'
    fid = fopen(cfile.name,'rt') ;
    X = fread(fid) ;
    fclose(fid) ;
    X = char(X.') ;
    % replace string S1 with string S2
    Y = strrep(X, 'double p,', 'double *p,') ;
    fid2 = fopen(cfile.name,'wt') ;
    fwrite(fid2,Y) ;
    fclose (fid2) ;
end

% mex ... ;
if choiceMani(1)
    mex fold_Manifold.c
    mex fold_NV.c
end
if choiceMani(2)
    mex modfold_Manifold.c;
    mex modfold_NV.c;
end
if choiceMani(3)
    mex hopf_Manifold.c;
    mex hopf_NV.c;
end
if choiceMani(4)
    mex modhopf_Manifold.c
    mex modhopf_NV.c
end
mex DDE_wrapped.c