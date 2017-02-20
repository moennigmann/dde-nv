%% Test Script to test the class EqualityConstraint


%clean up
clear all
close all
clc

% create variables to fill with
vars=VariableVector([2;3;4;5],0,[{'w'};{'x'};{'y'};{'z'}]);

% define rhs of equalities
conFun=@(a)[a(1)^2-1;a(2)*a(1)+1;a(3)];

%create instance

anEqCon=EqualityConstraint(conFun,3,vars,0);

% and remove old stuff again
clear conFun vars

% show indices prior to shifting
disp('eqIndex')
disp(anEqCon.eqIndex)
disp('varIndex')
disp(anEqCon.vars.index)

% shift up
anEqCon.shiftIndex(7,3);

disp('eqIndex')
disp(anEqCon.eqIndex)
disp('varIndex')
disp(anEqCon.vars.index)


% shift down
anEqCon.shiftIndex(-4,3);

disp('eqIndex')
disp(anEqCon.eqIndex)
disp('varIndex')
disp(anEqCon.vars.index)


%shift to much down to generate errors
try
    anEqCon.shiftIndex(-4,-3);
catch err
    disp(err.message)
end

try
    anEqCon.shiftIndex(-3,-13);
catch err
    disp(err.message)
end

