%% Demo of DDE-NV
% The following script demonstrates the basic steps to run a steady state
% optimization of a delayed system with robust exponential stability
% constraints using the DDE-NV toolbox


%% prepare matlab
close all
clear
clc

%% define system to begin with
a=5;
g=2;

param=[a;g];
xGue=[1;1];

% name states and parameters
xNom=VariableVector(xGue,0,[{'juvenile'};{'mature'}]);
alphaNom=VariableVector(param,2,[{'a'};{'g'}]);
p=VariableVector([],4,[]);

% define state dependent delay
tau=@(x,alpha,p)1.5-0.5*exp(-x(2));

% define system dynamics
popModel=@(x,xtau,alpha,p)[alpha(1)*x(2)-alpha(2)*x(1)-alpha(1)*exp(-alpha(2)*tau(x,alpha,p))*xtau(2,1);
    alpha(1)*exp(-alpha(2)*tau(x,alpha,p))*xtau(2,1)-x(2)^2];


% collect them in an object 
popDDE=DDE(@(x,xtau,alpha,p)popModel(x,xtau,alpha,p)',tau,xNom,alphaNom,p);


clear popModel
clear tau
clear alphaNom
clear xGue
clear param
clear a
clear g
clear gmaxStepLength

%% construct object ''aDDENLP''
% combine system dynamics, cost functions, box contraints etc. to form an
% nonlinear programm

J=@(x)x(1)/(x(1)+x(2));

aDDENLP=DDENLP(J,popDDE,xNom,[0;0],[Inf;Inf],[0;0],[10;4],-Inf(0,1),Inf(0,1));

clear popDDE
clear J
clear xNom

% define numerical options
aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-9,'TolX',1e-9);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-6,'RelLineSrchBnd',1e-5);

aDDENLP.maxAllowedRealPart=-0.1;

%% initialize steady state constraints
aDDENLP.initializeStSt();

%% add normal vector constraints
% these are responsible for the robust exponential stability

% compile c-code describing the constraints
mex populationModelModFoldMani.c
mex populationModelModFoldNV.c

% normal vector constraint for first exponential stability boundary
alphaCrit1=VariableVector([5;0.1],0,[{'acrit1'};{'gcrit1'}]);
xCrit1=VariableVector([29.87;4.3],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1


% normal vector constraint for second exponential stability boundary
alphaCrit2=VariableVector([4.6287;3.7225],0,[{'acrit1'};{'gcrit1'}]);
xCrit2=VariableVector([0.1144;0.0939],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit2,alphaCrit2,p);
aDDENLP.NVCon(end).vars.w1.values=[0.7817; 0.6236];
clear alphaCrit2
clear xCrit2

%% initatialize constraints and prepare optimization
aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();

%% visualize the  manifolds

manifoldSlices(1)=ManifoldSlice(aDDENLP.NVCon(1),[1 2]);
manifoldSlices(1).initStepLength=0.0005;
manifoldSlices(1).maxStepLength(1) = 0.05;
manifoldSlices(2)=ManifoldSlice(aDDENLP.NVCon(2),[1 2]);
manifoldSlices(2).initStepLength=0.0005;
manifoldSlices(2).maxStepLength(1) = 0.05;

figure(1);
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'ko');



hold on
box on
grid on
axis equal
xlim([0 8])
ylim([0 5])
[x,y] = circle(aDDENLP.minDist*sqrt(aDDENLP.nAlpha),aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2));
plot(x,y,'k');

manifoldSlices.maniContin2DbothDirections(150);

plot(aDDENLP.vars.critical(1).alpha.values(1),aDDENLP.vars.critical(1).alpha.values(2),'kx');
plot(aDDENLP.vars.critical(2).alpha.values(1),aDDENLP.vars.critical(2).alpha.values(2),'kx');
plot(manifoldSlices);


%% run optimization

aDDENLP.runOptim();
aDDENLP.deconstructOptimum();

% plot optimal nominal point
plot(aDDENLP.vars.critical(1).alpha.values(1),aDDENLP.vars.critical(1).alpha.values(2),'kx');
plot(aDDENLP.vars.critical(2).alpha.values(1),aDDENLP.vars.critical(2).alpha.values(2),'kx');
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'ko');
[x,y] = circle(aDDENLP.minDist*sqrt(aDDENLP.nAlpha),aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2));
plot(x,y,'k');

%% run simulation
% simulation of system dynamics with optimal parameters

point=aDDENLP.vars.nominal(end);

history=abs(point.x.values).*[0.9;1.1];
options=ddeset();
sol=ddesd(aDDENLP,point,history,[0 30],options);

figure(2);clf;
plot(sol.x,[sol.y(1,:)/point.x.values(1);sol.y(2,:)/point.x.values(2)])

%% check stability
% calculate eigenvalues at special points of interest

[~]=aDDENLP.checkStabilityPoint('nominal');
[~]=aDDENLP.checkStabilityPoint('critical');

[~]=aDDENLP.checkStabilityAtVertices();

[~]=aDDENLP.checkStabilityAtRandom();
