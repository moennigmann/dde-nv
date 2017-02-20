%% Test Script to test the class VariableVector


%clean up
clear all
close all
clc

%% test naming

disp('NAMES TEST');

try
    vars=VariableVector([2;3;4;5],0);
catch err
    disp(err.message)
end

disp(vars.names)


try
    vars=VariableVector([2;3;4;5],0,{'bla'});
catch err
    disp(err.message)
end

disp(vars.names)


clear vars

%% test the shifting method
disp('SHIFT TEST')


vars=VariableVector([2;3;4;5],1,[{'w'};{'x'};{'y'};{'z'}]);

disp(vars.index)
vars.shiftIndex(-1)
disp(vars.index)
try
    vars.shiftIndex(-1)
    disp(vars.index)
catch err
    disp(err.message)
end

