
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

%% define system dynamics
tau=@(x,p)1.5-0.5*exp(-x(2));
popModel=@(x,xtau,p)[p(1)*x(2)-p(2)*x(1)-p(1)*exp(-p(2)*tau(x,p))*xtau(2,1);
    p(1)*exp(-p(2)*tau(x,p))*xtau(2,1)-x(2)^2];

%% collect them in an object 
popDDE=DDE(@(x,xtau,p)popModel(x,xtau,p)',tau,xNom,alphaNom);


clear popModel
clear tau
clear alphaNom
clear xGue
clear param
clear a
clear g

%% construct object ''aDDENLP''

J=@(x)10*x(1)/(x(1)+x(2));

aDDENLP=DDENLP(J,popDDE,xNom,[0;0],[Inf;Inf],[0;0],[10;4]);
clear popDDE
clear J
clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-9);


%% initialize steady state constraints
aDDENLP=initializeStSt(aDDENLP);

%% add NV Cons

alphaCrit1=VariableVector([5;0.1],0,[{'acrit1'};{'gcrit1'}]);
xCrit1=VariableVector([29.87;4.3],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit1,alphaCrit1);
clear alphaCrit1
clear xCrit1

alphaCrit2=VariableVector([4.6287;3.7225],0,[{'acrit1'};{'gcrit1'}]);
xCrit2=VariableVector([0.1144;0.0939],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit2,alphaCrit2);
aDDENLP.NVCon(2).vars.w1.values=[0.7817; 0.6236];
clear alphaCrit2
clear xCrit2


%% initatialize constraints and prepare optimization
aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();

aDDENLP.fixedUncertParamIndex=[];

%% run optimization

aDDENLP.runOptim();
aDDENLP.deconstructOptimum();

%% run simulation

point=aDDENLP.vars.critical(2);

history=abs(point.x.values).*[0.9;1.1];
options=ddeset();
sol=ddesd(aDDENLP,point,history,[0 30],options);

figure(2);clf;
% plot(sol.x,sol.y);
plot(sol.x,[sol.y(1,:)/point.x.values(1);sol.y(2,:)/point.x.values(2)])

%% check stability

[~,aDDENLP]=checkStabilityPoint(aDDENLP,'nominal');
[~,aDDENLP]=checkStabilityPoint(aDDENLP,'critical');

[~,aDDENLP]=checkStabilityAtVertices(aDDENLP);

[~,aDDENLP]=checkStabilityAtRandom(aDDENLP);

