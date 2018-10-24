%% clear all
close all
clear
clc

%% name states and parameters
xNom = VariableVector([ 10,0,10,0,10,0]',0,[{'p1'};{'q1'};{'p2'};{'q2'};{'p3'};{'q3'}]);
p = VariableVector([],0,[]);
% This point is already unstable:
alphaNom = VariableVector([2.9 2.9 3.3 3.1 4 2.9 60]',6,[{'hD1'};{'hD2'};{'hD3'};{'hP1'};{'hP2'};{'hP3'};{'d'};]);



%% Define System dynamics

mySCmodel = @supplyChainModel;


%% collect them in an object

scDDE=DDE(@(x,xtau,alpha,pvv)mySCmodel(x,xtau,alpha,pvv)',@supplyChainDelays,xNom,alphaNom,p);

scDDE.hopfManiHandle = @supplyChainHopfMani;
scDDE.hopfNVHandle = @supplyChainHopfNV;

%% construct object ''aDDENLP''

% Cost function
J=@(x)(-x(6+7));

% set boundaries
ubx = [Inf;Inf;Inf;Inf;Inf;Inf];
uba = [Inf;Inf;Inf;Inf;Inf;Inf;1e3];
lbx = [-Inf;-Inf;-Inf;-Inf;-Inf;-Inf];
lba = [2.8; 2.8; 3.2; 2.75; 3.9; 2.65; sqrt(7)];

aDDENLP=DDENLP(J,scDDE,xNom,lbx,ubx,lba,uba,-Inf(0,1),Inf(0,1));

clear scDDE
clear J
clear xNom lba uba lbx ubx

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',1000000,'display','off','TolFun',1e-18,'TolX',1e-12);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-12,'TolCon',1e-15);
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',100,'MaxFunEvals',200000,'display','iter','TolFun',1e-12,'TolCon',1e-12,'RelLineSrchBnd',1);


aDDENLP.allowedEigsInClosedRightHP=2;

%% initialize steady state constraints
aDDENLP.initializeStSt();
aDDENLP.checkStabilityPoint('nominal');

%% prepare optimization
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();

%% run optimization
aDDENLP.runOptimAddingNewManifolds(3);
aDDENLP.deconstructOptimum();

aDDENLP.checkStabilityPoint('nominal');

%% run simulation

options = ddeset('AbsTol',1e-6,'RelTol',1e-4);
point=aDDENLP.vars.nominal(1);

history=abs(point.x.values).*[0.999;1.111;1;1;1;1];
sol=aDDENLP.ddesd(aDDENLP.vars.nominal, history, [0,1000],options);

figure(2);clf;
plot(sol.x,sol.y(5,:))