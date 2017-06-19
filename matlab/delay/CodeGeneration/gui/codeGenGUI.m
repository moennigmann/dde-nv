%% Program to open input dialogs (gui), code generation, code manipulation and compiling
% add folders gui, explicitAEs and codegen to matlab

%% start gui

systemInput;
% systemInput_limit;

% if additional equations are needed (eg:system equations have additional
% parameter)
explicitAEs;

%% generate code for maple

codeGeneration;
%%
clear aes alphavec ans answer anum del delnum fileID i j l tauname way x xdot xnames xnum

system('maple codegen.txt');

%% code manipulation to make it compilable for matlab

if y(1) % Fold
    mexWrapperGen_Fold;
    mexWrapperGen_FoldNV;
end
if y(2) % modFold
    mexWrapperGen_Modfold;
    mexWrapperGen_ModfoldNV;
end
if y(3) % Hopf
    mexWrapperGen_HopfMani;
    mexWrapperGen_HopfNV;
end
if y(4) % modHopf
    mexWrapperGen_ModHopf;
    mexWrapperGen_ModHopfNV;
end

%% Compile for matlab

% mex ... ;
if y(1)
    mex fold_Manifold.c
    mex fold_NV.c
end
if y(2)
    mex modfold_Manifold.c;
    mex modfold_NV.c;
end
if y(3)
    mex hopf_Manifold.c;
    mex hopf_NV.c;
end
if y(4)
    mex modhopf_Manifold.c
    mex modhopf_NV.c
end