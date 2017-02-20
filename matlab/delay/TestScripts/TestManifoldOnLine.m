
%% prepare matlab
close all
clear all
clc

%% define system to begin with
a=5;
g=2;

param=[a;g];
xGue=[1;1];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'juvenile'};{'mature'}]);
alphaNom=VariableVector(param,2,[{'a'};{'g'}]);
p=VariableVector([-0.0;-0.7],4,[{'bla'};{'blabla'}]);

%% define system dynamics
tau=@(x,alpha,p)1.5-0.5*exp(-x(2));
popModel=@(x,xtau,alpha,p)[alpha(1)*x(2)-alpha(2)*x(1)-alpha(1)*exp(-alpha(2)*tau(x,alpha,p))*xtau(2,1);
    alpha(1)*exp(-alpha(2)*tau(x,alpha,p))*xtau(2,1)-x(2)^2];


%% collect them in an object 
popDDE=DDE(@(x,xtau,alpha,p)popModel(x,xtau,alpha,p)',tau,xNom,alphaNom,p);


clear popModel
clear tau
clear alphaNom
clear xGue
clear param
clear a
clear g

%% construct object ''aDDENLP''

J=@(x)x(1)/(x(1)+x(2));

aDDENLP=DDENLP(J,popDDE,xNom,[0;0],[Inf;Inf],[0;0],[10;4],[-3.5;-4],[Inf;Inf]);
clear popDDE
clear J
clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-9);


%% initialize steady state constraints
aDDENLP.initializeStSt();

%% add NV Cons

% alphaCrit1=VariableVector([5;0.1],0,[{'acrit1'};{'gcrit1'}]);
% xCrit1=VariableVector([29.87;4.3],0,{'xcrit1'});
% 
% aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit1,alphaCrit1,p);
% clear alphaCrit1
% clear xCrit1
% 
% alphaCrit2=VariableVector([4.6287;3.7225],0,[{'acrit1'};{'gcrit1'}]);
% xCrit2=VariableVector([0.1144;0.0939],0,{'xcrit1'});
% 
% aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit2,alphaCrit2,p);
% aDDENLP.NVCon(end).vars.w1.values=[0.7817; 0.6236];
% clear alphaCrit2
% clear xCrit2
% 
point1=aDDENLP.vars.nominal;
point2.x=VariableVector([0.5,0.5]',Inf,[{'juv'},{'mat'}]');
point2.alpha=VariableVector([2,4]',Inf,[{'a'},{'g'}]');
point2.p=VariableVector([-0.0;-0.7],Inf,[{'bla'};{'blabla'}]);
intermediatePoint = findManifoldPointOnLine(aDDENLP,'modfold',@populationModelModFoldMani,point1,point2);

