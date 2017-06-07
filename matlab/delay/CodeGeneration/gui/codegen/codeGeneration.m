%% Code generation for maple - control
% to do: codeGenModHopf, codeGenFold
% to do: mexWrapper_ModHopdf _ModHopfNV _Fold _FoldNV


% initialisation
codeGenInit;
% get p vector (later not necessary)
% getp;
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