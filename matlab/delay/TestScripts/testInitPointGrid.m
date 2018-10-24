%% Optimization for multiple initial points
% a grid of initial points is created to notice convergence to local
% optima


%% clean up
close all
clear
clc



%% prepare matlab
close all
clear
clc

%% define system to begin with
a=5;
g=2;

param=[a;g];

% multiParam=param+2*randn(2,20)


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
% clear param
clear a
clear g

%% construct object ''aDDENLP''

J=@(x)sin(x(3)).*cos(x(4));


[agrid,ggrid]=meshgrid(linspace(-1,11,100),linspace(-1,5,100));

Jgrid=NaN(size(agrid));
for jj=1:numel(agrid)
    Jgrid(jj)=J([0,0,agrid(jj),ggrid(jj)]);
end


aDDENLP=DDENLP(J,popDDE,xNom,[0;0],[Inf;Inf],[0;0],[10;4],-Inf(0,1),Inf(0,1));
clear popDDE
clear J
clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-9,'TolX',1e-9);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-6,'RelLineSrchBnd',1e-5);




%% run optimization
aDDENLP.maxAllowedRealPart=-0.1;


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

%% specific stuff for different initial points

figure(1);
hold on;
axis equal;
box on;
grid on;
xlim([0 20]);

aDDENLP.runOptimMultipleInitPoints([param-1, param, param+1]);



figure(1)

plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'.')

[uncertRegionAlpha, uncertRegionGamma] = circle( sqrt(2), aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2) );

plot(uncertRegionAlpha, uncertRegionGamma);

mani=load('popCritMani1');
plot(mani.alpha,mani.gamma,[0 20],[0.1,0.1])
contour(agrid,ggrid,Jgrid)


[~]=aDDENLP.checkStabilityPoint('nominal');
[~]=aDDENLP.checkStabilityPoint('critical');
[~]=aDDENLP.checkStabilityAtVertices()
[~]=aDDENLP.checkStabilityAtRandom();

