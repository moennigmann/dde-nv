
%% prepare matlab
close all
clear all
clc


a=5;
g=2;

param=[a;g];
xGue=[1;1];
p=[1];

xNom=VariableVector(xGue,0,[{'juvenile'};{'mature'}]);
alphaNom=VariableVector(param,2,[{'a'};{'g'}]);
pNom=VariableVector(p,4,[{'bla'}]);

tau=@(x,p)1.5-0.5*exp(-x(2));
popModel=@(x,xtau,alpha,p)[alpha(1)*x(2)-alpha(2)*x(1)-alpha(1)*exp(-alpha(2)*tau(x,alpha))*xtau(2,1);
    alpha(1)*exp(-alpha(2)*tau(x,alpha))*xtau(2,1)-x(2)^2];


popDDE=DDE(@(x,xtau,alpha,p)popModel(x,xtau,alpha,p)',tau,xNom,alphaNom,pNom);
nvVars1=varCollection('modfold',4+length(p),xNom,alphaNom,pNom);
nvVars2=varCollection('modhopf',18+length(p),xNom,alphaNom,pNom);


nvVars1.x.values = [0.188230344817496; 0.0916303916372130];
nvVars1.alpha.values = [9.13026377952101; 4.40000000000000];
nvVars1.w1.values = [1; 0];

nvVars2.alpha.values=[5;0.1];

clear popModel
clear tau
clear xGue
clear param
clear a
clear g

%% try DDENLP

J=@(x)x(1)/(x(1)+x(2));

aDDENLP=DDENLP(J,popDDE,xNom,[0.0;0.0],[Inf;Inf],[4;2],[10;3]);
clear popDDE
clear J
clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','display','off');


aDDENLP=initializeStSt(aDDENLP);



%% construct a normal vector constraint



myNVCon1 = NVConstraint(aDDENLP,'modfold', @populationModelModFoldMani, @populationModelModFoldNV, nvVars1);
myNVCon1 = findManifoldPoint(myNVCon1,nvVars1);
myNVCon1 = findClosestCriticalPoint(myNVCon1,alphaNom);
myNVCon1 = findNormalVector(myNVCon1);
myNVCon1 = findConnection(myNVCon1,alphaNom);

myNVCon2 = NVConstraint(aDDENLP,'modfold', @populationModelModFoldMani, @populationModelModFoldNV, nvVars2);
myNVCon2 = findManifoldPoint(myNVCon2, nvVars2);
myNVCon2 = findClosestCriticalPoint(myNVCon2,alphaNom);
myNVCon2 = findNormalVector(myNVCon2);
myNVCon2 = findConnection(myNVCon2,alphaNom);

