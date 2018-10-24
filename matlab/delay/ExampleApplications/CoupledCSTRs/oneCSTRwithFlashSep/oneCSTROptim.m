%% Optimization of CSTR with Flash Separator and Recycling and rection A->B->C
 
%% prepare matlab
close all
clear all

% clear screen only if this script is called directly and skip screen
% clearance when script was called by another script
functionStackVar=dbstack;

if strcmp(mfilename,functionStackVar(end).name)
    clc
end
clear functionStackVar

%% define system to begin with
Fr=50; % m^3/h

T01=400; %K
F01=60;% m^3/h
Tc1=480; %K
V1=1*1000; %dm^3

TcF=480; %K
VF=1*1000; %dm^3

beta=4.5*10;

param=[Fr, beta, T01, F01, Tc1, V1, TcF, VF]';
xGue=[0.433749162609543;0.514872752267692;5.215504372524927e+02;0.291959548082289;0.632458648491924;4.802802298036649e+02];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'cA1'};{'cB1'};{'T1'};{'cAF'};{'cBF'};{'TF'}]);
alphaNom=VariableVector(param,6,[{'Fr'}; {'beta'};{'T01'};{'F01'};{'Tc1'};{'V1'};{'TcF'};{'VF'}]);
p=VariableVector([],0,[]);


clear Fr T01 F01 Tc1 V1 TcF VF beta xGue param


%% define system dynamics
tau=@(x,alpha,p)[5*350/(alpha(1));
    5*350/(alpha(1)+alpha(4))];

cstrModel=@(x,xtau,alpha,p)0.01*oneCSTRoneFlashSepDDE(x,xtau,alpha,p);



%% collect them in an object 
cstrDDE=DDE(@(x,xtau,alpha,p)cstrModel(x,xtau,alpha,p),tau,xNom,alphaNom,p);


clear cstrModel
clear tau
clear alphaNom
%% construct object ''aDDENLP''

% J=@costFun1;
J=@(x)100*returnOnInvest1(x);

setNCSTRBoundaries

oneCSTRDDENLP=DDENLP(J,cstrDDE,xNom,...
    [Cmin;Cmin;Tmin;Cmin;Cmin;Tmin],...
    [Cmax;Cmax;Tmax;Cmax;Cmax;Tmax],...
    [Frmin, betamin, Tmin, F0min, Tcmin, Vmin, Tcmin, Vmin]'+minDist*sqrt(8),...
    [Frmax, betamax, T0max, F0max, Tmax, Vmax, Tmax, Vmax]'-minDist*sqrt(8),...
    -Inf(0,1),Inf(0,1));

clear cstrDDE J xNom Cmin Tmin Cmax Tmax Frmin betamin F0min Vmin Frmax betamax Tmax F0max Vmax

oneCSTRDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off');
oneCSTRDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-6,'TolX',1e-4);
oneCSTRDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-9);


%% initialize steady state constraints
oneCSTRDDENLP.initializeStSt();

%% add NV Cons

alphaCrit1=VariableVector([53.871204633059900 45 400 65.398303761840860 480 1000 480 1000]',0,[{'alphacrit1'}]);
xCrit1=VariableVector([0.435362112802818;0.512881904267267;5.251569436184192e+02;0.294258135643240;0.629781199583491;4.803300298743896e+02],0,{'xcrit1'});

oneCSTRDDENLP.addNVCon('hopf',@oneCSTRoneFlashSepHopfMani,@oneCSTRoneFlashSepHopfManiNV,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1

oneCSTRDDENLP.NVCon.vars.omega.values=0.3214762433635;
oneCSTRDDENLP.NVCon.vars.w1.values=real([-1.655225591103498e-04 + 8.417850970913898e-04i;1.598586611723015e-04 - 7.284243731761511e-04i;0.999972774717766 - 0.000000082788724i;2.857399109998191e-05 - 8.309278478180857e-05i;-2.672920252770429e-05 + 7.160473420218546e-05i;-0.000314016144977 - 0.007283234021262i]);
oneCSTRDDENLP.NVCon.vars.w2.values=imag([-1.655225591103498e-04 + 8.417850970913898e-04i;1.598586611723015e-04 - 7.284243731761511e-04i;0.999972774717766 - 0.000000082788724i;2.857399109998191e-05 - 8.309278478180857e-05i;-2.672920252770429e-05 + 7.160473420218546e-05i;-0.000314016144977 - 0.007283234021262i]);


%% initatialize constraints and prepare optimization
oneCSTRDDENLP.initNVCons();

oneCSTRDDENLP.concatConstraints();
oneCSTRDDENLP.concatInitPoints();

oneCSTRDDENLP.fixedUncertParamIndex=[];

%% run optimization

oneCSTRDDENLP.runOptim();
oneCSTRDDENLP.deconstructOptimum();

%% run simulation

point=oneCSTRDDENLP.vars.nominal(1);

history=abs(point.x.values).*[0.9;1.1;1;1;1;1];
options=ddeset();
sol=ddesd(oneCSTRDDENLP,point,history,[0 10000],options);

figure(2);clf;
hold on
plot(sol.x,sol.y(1:3:end,:),'-.')
plot(sol.x,sol.y(2:3:end,:))
plot([sol.x(1,1),sol.x(1,end)],[point.x.values(2:3:end),point.x.values(2:3:end)]','--')

%% check stability

oneCSTRDDENLP.numMinEig=-0.5;

maxEigN=oneCSTRDDENLP.checkStabilityPoint('nominal');
maxEigC=oneCSTRDDENLP.checkStabilityPoint('critical');

% [~,oneCSTRDDENLP]=checkStabilityAtVertices(oneCSTRDDENLP);
% 
% [~,oneCSTRDDENLP]=checkStabilityAtRandom(oneCSTRDDENLP);

