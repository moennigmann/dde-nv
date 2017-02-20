
%% clean up workspace
close all
clear all
clc

% define parameters
Fr=50; % m^3/h
beta=4.5*10;

T01=400; %K
F01=60;% m^3/h
Tc1=480; %K
V1=1*1000; %dm^3

TcF=480; %K
VF=1*1000; %dm^3

param=[Fr,beta,...
    T01,F01,Tc1,V1,...
    TcF,VF]';

xGue=[0.433749162609543;0.514872752267692;5.215504372524927e+02;0.291959548082289;0.632458648491924;4.802802298036649e+02];

xNom=VariableVector(xGue,0,[{'cA1'},{'cB1'},{'T1'},{'cA2'},{'cB2'},{'T2'}]');
alphaNom=VariableVector(param,6,[{'Fr'},{'beta'},{'T01'},{'F01'},{'Tc1'},{'V1'},...
    {'TcF'},{'VF'}]');

aDDE=DDE(@(x,xtau,p)oneCSTRoneFlashSepDDE(x,xtau,p)',@(x,p)NCSTRDelayVector(x,p,2),xNom,alphaNom);

nomVars.x=xNom;
nomVars.alpha=alphaNom;

aStStCon=StStConstraint(aDDE,nomVars);

res=aStStCon.conFun([xNom.values;alphaNom.values]);