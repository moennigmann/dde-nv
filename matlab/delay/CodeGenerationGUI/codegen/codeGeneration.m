%> @file codeGeneration.m
%> @brief main code for generating the maple code. first the name of the
%> txt file is defined with the name of the project


%% Code generation for maple - control

% define name for code generation files
cg_name = ['mapleCodeGeneration_' p_name '.txt'];

% initialisation
codeGenInit;

% define system
codeGenDefSys;

% pick relevant equations of normal vector system - Fold
if choiceMani(1)
    codeGenFold;
end
% pick relevant equations of normal vector system - Manifold modfold
if choiceMani(2)
    codeGenModFold;
end
% pick relevant equations of normal vector system - Hopf
if choiceMani(3)
    codeGenHopf;
end
if choiceMani(4)
    codeGenModHopf;
end
% get dde equations
codeGenDDE;