
%% prepare matlab
close all
clear
clc

%% define system to begin with
a=5;
g=2;

param=[a;g];
xGue=[1;1];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'juvenile'};{'mature'}]);
alphaNom=VariableVector(param,2,[{'a'};{'g'}]);
p=VariableVector([],4,[]);

%% define system dynamics
tau=@(x,alpha,p)1.5-0.5*exp(-x(2));
popModel=@(x,xtau,alpha,p)[alpha(1)*x(2)-alpha(2)*x(1)-alpha(1)*exp(-alpha(2)*tau(x,alpha,p))*xtau(2,1);
    alpha(1)*exp(-alpha(2)*tau(x,alpha,p))*xtau(2,1)-x(2)^2];


%% collect them in an object 
popDDE=DDE(@(x,xtau,alpha,p)popModel(x,xtau,alpha,p)',tau,xNom,alphaNom,p);

popDDE.modfoldManiHandle = @populationModelModFoldMani;
popDDE.modfoldNVHandle = @populationModelModFoldNV;


clear popModel
clear tau
clear alphaNom
clear xGue
clear param
clear a
clear g

%% construct object ''aDDENLP''



J=@(x)x(1)/(x(1)+x(2));

aDDENLP=DDENLP(J,popDDE,xNom,[0.05;0.05],[Inf;Inf],[0.01;-5],[10;4],[],[]);
clear popDDE
clear J
clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',500,'MaxFunEvals',200000,'display','final','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-6,'RelLineSrchBnd',1e-5);


%% initialize steady state constraints
aDDENLP.initializeStSt();

%% add NV Cons

alphaCrit1=VariableVector([5;0.1],0,[{'acrit1'};{'gcrit1'}]);
xCrit1=VariableVector([29.87;4.3],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1
 
% alphaCrit2=VariableVector([4.6287;3.7225],0,[{'acrit1'};{'gcrit1'}]);
% xCrit2=VariableVector([0.1144;0.0939],0,{'xcrit1'});
% 
% aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit2,alphaCrit2,p);
% aDDENLP.NVCon(end).vars.w1.values=[0.7817; 0.6236];
% clear alphaCrit2
% clear xCrit2


figure(1);
hold on;
axis equal;
box on;
grid on;
xlim([0 20]);

mani=load('popCritMani1');
plot(mani.alpha,mani.gamma,[0 20],[0.1,0.1])

for ii = 1:length(aDDENLP.vars.critical)
    alphaCrit = aDDENLP.vars.critical(ii).alpha.values(1);
    gammaCrit = aDDENLP.vars.critical(ii).alpha.values(2);
    plot(alphaCrit,gammaCrit,'x')
end


plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'.')


%% initatialize constraints and prepare optimization
aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();


%% run optimization
aDDENLP.maxAllowedRealPart=-0.1;

aDDENLP.runOptimAddingNewManifolds(10);


figure(1);
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'o')
plot(aDDENLP.vars.critical(end).alpha.values(1),aDDENLP.vars.critical(end).alpha.values(2),'>')

[uncertRegionAlpha, uncertRegionGamma] = circle( sqrt(2), aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2) );

plot(uncertRegionAlpha, uncertRegionGamma);
 
%% run simulation

point=aDDENLP.vars.critical(end);

history=abs(point.x.values).*[0.9;1.1];
options=ddeset();
sol=ddesd(aDDENLP,point,history,[0 30],options);

figure(2);clf;
plot(sol.x,[sol.y(1,:)/point.x.values(1);sol.y(2,:)/point.x.values(2)])

%% check stability

[~]=aDDENLP.checkStabilityPoint('nominal');
[~]=aDDENLP.checkStabilityPoint('critical');

[~]=aDDENLP.checkStabilityAtVertices();

[~]=aDDENLP.checkStabilityAtRandom();
