
%% prepare matlab
close all
clear
clc

%% define system to begin with
a=0.5;
g=0.2;

a =5;
g=4;
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

aDDENLP=DDENLP(J,popDDE,xNom,[-Inf;-Inf],[Inf;Inf],[0.01;-5],[40;10],[],[]);
clear popDDE
clear J
clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',500,'MaxFunEvals',200000,'display','iter','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-6,'RelLineSrchBnd',1e-5);


%% initialize steady state constraints
aDDENLP.initializeStSt();
[~]=aDDENLP.checkStabilityPoint('nominal');


%% add NV Cons

alphaCrit1=VariableVector([5;0.1],0,[{'acrit1'};{'gcrit1'}]);
xCrit1=VariableVector([29.87;4.3],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1
 
alphaCrit2=VariableVector([4.6287;3.7225],0,[{'acrit1'};{'gcrit1'}]);
xCrit2=VariableVector([0.1144;0.0939],0,{'xcrit1'});

aDDENLP.addNVCon('modfold',@populationModelModFoldMani,@populationModelModFoldNV,xCrit2,alphaCrit2,p);
aDDENLP.NVCon(end).vars.w1.values=[0.7817; 0.6236];
clear alphaCrit2
clear xCrit2


figure(1);
hold on;
axis equal;
box on;
grid on;
xlim([0 20]);

mani=load('popCritMani1');
plot(mani.alpha,mani.gamma,[0 20],[0.1,0.1])

plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'.')

[uncertRegionAlpha, uncertRegionGamma] = circle( sqrt(2), aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2) );

plot(uncertRegionAlpha, uncertRegionGamma);


%% initatialize constraints and prepare optimization
aDDENLP.initNVCons();

plot(aDDENLP.vars.critical(2).alpha.values(1),aDDENLP.vars.critical(2).alpha.values(2),'x')
quiver(aDDENLP.vars.critical(2).alpha.values(1),aDDENLP.vars.critical(2).alpha.values(2),aDDENLP.vars.critical(2).r.values(1),aDDENLP.vars.critical(2).r.values(2))
plot(aDDENLP.vars.critical(1).alpha.values(1),aDDENLP.vars.critical(1).alpha.values(2),'x')
quiver(aDDENLP.vars.critical(1).alpha.values(1),aDDENLP.vars.critical(1).alpha.values(2),aDDENLP.vars.critical(1).r.values(1),aDDENLP.vars.critical(1).r.values(2))


aDDENLP.moveAwayFromManifolds(0.5, 2, 100);
return
fprintf('alphaNom = \n [%s]''\n', strjoin(cellstr(num2str(aDDENLP.vars.nominal.alpha.values(:))),', '));
plot(aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2),'o')

plot(aDDENLP.vars.critical(1).alpha.values(1), aDDENLP.vars.critical(1).alpha.values(2),'>')
plot(aDDENLP.vars.critical(2).alpha.values(1), aDDENLP.vars.critical(2).alpha.values(2),'>')

[uncertRegionAlpha, uncertRegionGamma] = circle( sqrt(2), aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2) );

plot(uncertRegionAlpha, uncertRegionGamma);


aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();


%% run optimization
aDDENLP.maxAllowedRealPart=-0.1;

aDDENLP.runOptim();
aDDENLP.deconstructOptimum;


figure(1);
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'o')
plot(aDDENLP.vars.critical(1).alpha.values(1),aDDENLP.vars.critical(1).alpha.values(2),'>')
plot(aDDENLP.vars.critical(end).alpha.values(1),aDDENLP.vars.critical(end).alpha.values(2),'>')

[uncertRegionAlpha, uncertRegionGamma] = circle( sqrt(2), aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2) );

plot(uncertRegionAlpha, uncertRegionGamma);
 
% %% run simulation
% 
% point=aDDENLP.vars.critical(end);
% 
% history=abs(point.x.values).*[0.9;1.1];
% options=ddeset();
% sol=ddesd(aDDENLP,point,history,[0 30],options);
% 
% figure(2);clf;
% plot(sol.x,[sol.y(1,:)/point.x.values(1);sol.y(2,:)/point.x.values(2)])
% 
% %% check stability
% 
% [~]=aDDENLP.checkStabilityPoint('nominal');
% [~]=aDDENLP.checkStabilityPoint('critical');
% 
% [~]=aDDENLP.checkStabilityAtVertices();
% 
% [~]=aDDENLP.checkStabilityAtRandom();
