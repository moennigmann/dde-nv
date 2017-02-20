%% Optimization of CSTR with Flash Separator and Recycling and rection A->B->C
 
%% prepare matlab
close all
clear all
clc

%% define system to begin with
Fr=5; % m^3/h
beta=5*10;

T01=400; %K
F01=5;% m^3/h
Tc1=480; %K
V1=2*1000; %dm^3

T02=405; %K
F02=5;% m^3/h
Tc2=485; %K
V2=1.5*1000; %dm^3


T03=410; %K
F03=5;% m^3/h
Tc3=490; %K
V3=1*1000; %dm^3


TcF=480; %K
VF=1*1000; %dm^3

param1=[Fr, beta, T01, F01, Tc1, V1, T02, F02, Tc2, V2, T03, F03, Tc3, V3, TcF, VF]';
xGue1=[0.596646991868868;0.385230034504465;492.122444638444;0.495863392383575;0.477573246208705;500.122330492971;0.422755558869387;0.542125795048906;498.581136968797;0.294503618477956;0.655686307059929;480.067881934534];


param2=[3.46410521956844;86.9042406459036;293.464101615138;96.5358983848623;536.577139600221;103.464101604451;293.464101615138;96.5358983851387;614.625871007046;1668.98361145438;293.464101483153;3.46410161513775;653.799312251681;103.464104690120;338.464105564604;103.464101615138];
xGue2=[0.332854611445935;0.495598812890548;700;0.0394928150130609;0.246264591525011;699.999999212163;0.0328122368013648;0.239374694603774;676.480481170964;0.0307248160662403;0.237701417330125;347.992354235936];

lambda=0;

% lambda=0.718;% Hopf!


param=param1+lambda*(param2-param1);

% [1.80809086725597;58.8099444881404;323.552219308018;76.9045525499781;531.562810164766;638.209760190757;324.916650230040;70.3202199487999;577.887202445886;1621.34921710526;326.326679791659;2.70556574322655;605.535569558336;356.285246004059;380.025870926792;356.289256448028]

xGue=xGue1+lambda*(xGue2-xGue1);


%% name states and parameters
xNom=VariableVector(xGue,0,[{'cA1'};{'cB1'};{'T1'};{'cA2'};{'cB2'};{'T2'};{'cA3'};{'cB3'};{'T3'};{'cAF'};{'cBF'};{'TF'}]);
alphaNom=VariableVector(param,12,[{'Fr'}; {'beta'};{'T01'};{'F01'};{'Tc1'};{'V1'};{'T02'};{'F02'};{'Tc2'};{'V2'};{'T03'};{'F03'};{'Tc3'};{'V3'};{'TcF'};{'VF'}]);
p=VariableVector([],12+16,[]);



clear Fr T01 F01 Tc1 V1 T02 F02 Tc2 V2 T03 F03 Tc3 V3 TcF VF beta xGue param


%% define system dynamics
tau=@(x,alpha,p)[5*350/(alpha(1));
    5*350/(alpha(1)+alpha(4));
    5*350/(alpha(1)+alpha(4)+alpha(8))
    5*350/(alpha(1)+alpha(4)+alpha(8)+alpha(12))];

cstrModel=@(x,xtau,alpha,p)threeCSTRoneFlashSepDDE(x,xtau,alpha,p);



%% collect them in an object 
cstrDDE=DDE(@(x,xtau,alpha,p)0.01*cstrModel(x,xtau,alpha,p),tau,xNom,alphaNom,p);


clear cstrModel
clear tau
clear alphaNom


%% construct object ''aDDENLP''

J=@(x)100*returnOnInvest3(x);

setNCSTRBoundaries


threeCSTRDDENLP=DDENLP(J,cstrDDE,xNom,...
    [Cmin;Cmin;Tmin;Cmin;Cmin;Tmin;Cmin;Cmin;Tmin;Cmin;Cmin;Tmin],...
    [Cmax;Cmax;Tmax;Cmax;Cmax;Tmax;Cmax;Cmax;Tmax;Cmax;Cmax;Tmax],...
    [Frmin, betamin, Tmin, F0min, Tcmin, Vmin, Tmin, F0min, Tcmin, Vmin, Tmin, F0min, Tcmin, Vmin, Tcmin, Vmin]'+minDist*sqrt(12),...
    [Frmax, betamax, T0max, F0max, Tmax, Vmax, T0max, F0max, Tmax, Vmax, T0max, F0max, Tmax, Vmax, Tmax, Vmax]'-minDist*sqrt(12),...
    -Inf(0,1),...
    Inf(0,1));

clear cstrDDE J xNom Cmin Tmin Cmax Tmax Frmin betamin F0min Vmin Frmax betamax Tmax F0max Vmax



% threeCSTRDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000000,'MaxFunEvals',200000000,'display','off','TolFun',1e-7,'TolX',1e-7);
% threeCSTRDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',1000000,'MaxFunEvals',20000000,'display','iter-detailed','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9,'FunValCheck','on');
% threeCSTRDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',100000,'MaxFunEvals',2000000,'display','iter','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-3,'RelLineSrchBnd',0.02);

threeCSTRDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-5,'TolX',1e-7);
threeCSTRDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',20000,'MaxFunEvals',1000000,'display','iter','TolFun',1e-6,'TolX',1e-6,'TolCon',1e-6,'FunValCheck','on');
threeCSTRDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',100000,'MaxFunEvals',2000000,'display','iter','TolFun',1e-5,'TolX',1e-5,'TolCon',1e-1);






%% initialize steady state constraints
threeCSTRDDENLP.initializeStSt();


% 
% %% start a simulation
% 
% 
% figure(3);clf;
% point=threeCSTRDDENLP.vars.nominal(1);
% 
% history=abs(point.x.values).*[0.9;1.1;1;1;1;1;1;1;1;1;1;1];
% options=ddeset();
% sol=ddesd(threeCSTRDDENLP,point,history,[0 1000],options);
% 
% hold on
% plot(sol.x,sol.y(1:3:end,:))
% plot(sol.x,sol.y(2:3:end,:),'-.')
% plot([sol.x(1,1),sol.x(1,end)],[point.x.values(1:3:end),point.x.values(1:3:end)]','--')
% 
% 
% break


%% add NV Cons
% 
% alphaCrit1=VariableVector([5 8.25865471585134 400 2.49996668115040 480 1000 400 5 490 1000 400 5 500 1000 480 1000]',0,{'alphacrit1'});
% xCrit1=VariableVector([0.0409764337889367;0.335947788097982;501.682481215587;0.0385385954159356;0.301931065507839;532.903141158335;0.0292103579619715;0.252808868108509;547.792176122853;0.0129429175165232;0.218878166895381;480.396694512612],0,{'xcrit1'});
% 
% threeCSTRDDENLP.addNVCon('hopf',@threeCSTRoneFlashSepHopfMani,@threeCSTRoneFlashSepHopfManiNV,xCrit1,alphaCrit1,p);
% clear alphaCrit1
% clear xCrit1
% 
% threeCSTRDDENLP.NVCon.vars.omega.values=0.1283085055893;threeCSTRDDENLP.NVCon(2).vars.x.values
% threeCSTRDDENLP.NVCon.vars.w1.values=real([5.36950468855040e-08 - 4.10233886590683e-07i;
%     -9.35375999073203e-08 + 1.16972375406079e-07i;
%     -4.91756607293937e-05 - 0.000270361421746502i;
%     -1.02721518286624e-08 - 1.78578984562354e-08i;
%     -1.85568432493949e-09 + 5.36117702083144e-09i;
%     -7.19898588082646e-05 + 4.05654445591128e-05i;
%     -7.51946965900346e-05 + 0.000185110422543392i;
%     6.50589574401767e-05 + 4.18620409762054e-05i;
%     0.999983212874028 - 5.54206442237634e-05i;
%     7.11253395815784e-06 + 2.51046071912482e-06i;
%     1.41096752326775e-06 - 2.71003317970779e-06i;
%     0.00555516698133207 - 0.00160707972122930i]);
% threeCSTRDDENLP.NVCon.vars.w2.values=imag([5.36950468855040e-08 - 4.10233886590683e-07i;
%     -9.35375999073203e-08 + 1.16972375406079e-07i;
%     -4.91756607293937e-05 - 0.000270361421746502i;
%     -1.02721518286624e-08 - 1.78578984562354e-08i;
%     -1.85568432493949e-09 + 5.36117702083144e-09i;
%     -7.19898588082646e-05 + 4.05654445591128e-05i;
%     -7.51946965900346e-05 + 0.000185110422543392i;
%     6.50589574401767e-05 + 4.18620409762054e-05i;
%     0.999983212874028 - 5.54206442237634e-05i;
%     7.11253395815784e-06 + 2.51046071912482e-06i;
%     1.41096752326775e-06 - 2.71003317970779e-06i;
%     0.00555516698133207 - 0.00160707972122930i]);
% 



alphaCrit2=VariableVector(param1+0.69172*(param2-param1),0,{'alphacrit2'});
xCrit2=VariableVector(xGue1+0.69172*(xGue2-xGue1 ),0,{'xcrit2'});

threeCSTRDDENLP.addNVCon('hopf',@threeCSTRoneFlashSepHopfMani,@threeCSTRoneFlashSepHopfManiNV,xCrit2,alphaCrit2,p);
clear alphaCrit2
clear xCrit2

alphaCrit2=VariableVector(param1+0*0.7221*(param2-param1),0,{'alphacrit2'});
xCrit2=VariableVector(xGue1+0*0.7221*(xGue2-xGue1 ),0,{'xcrit2'});

xCrit2=VariableVector(threeCSTRDDENLP.vars.nominal.x.values,0,{'xcrit2'});

threeCSTRDDENLP.addNVCon('hopf',@threeCSTRoneFlashSepHopfMani,@threeCSTRoneFlashSepHopfManiNV,xCrit2,alphaCrit2,p);

% res=threeCSTRDDENLP.NVCon(2).eqAugSys.conFun([zeros(1,threeCSTRDDENLP.NVCon(2).vars.x.index(1)-1)...
% threeCSTRDDENLP.NVCon(2).vars.x.values',...
% threeCSTRDDENLP.NVCon(2).vars.alpha.values',...
% threeCSTRDDENLP.NVCon(2).vars.p.values',...
% 0,...
% threeCSTRDDENLP.NVCon(2).vars.x.values',...
% threeCSTRDDENLP.NVCon(2).vars.x.values']);
% res(1:12)


clear alphaCrit2
clear xCrit2


%% initatialize constraints and prepare optimization
% threeCSTRDDENLP.NVCon(1).findManifoldPoint(threeCSTRDDENLP.NVCon(1).vars)
bla1=threeCSTRDDENLP.NVCon(2).vars.x.values
threeCSTRDDENLP.NVCon(2).findManifoldPoint(threeCSTRDDENLP.NVCon(2).vars)
threeCSTRDDENLP.NVCon(2).findClosestCriticalPoint(threeCSTRDDENLP.vars.nominal.alpha)
bla2=threeCSTRDDENLP.NVCon(2).vars.x.values

bla1-bla2

threeCSTRDDENLP.numMinEig=-0.03;

% [maxEigN, eigsN] = threeCSTRDDENLP.checkStabilityPoint('nominal')
[maxEig, eigsC] = threeCSTRDDENLP.checkStabilityPoint('critical')

plot(eigsC,'*');

break
% threeCSTRDDENLP.NVCon.findNormalVector(threeCSTRDDENLP.vars.nominal.alpha);
% % threeCSTRDDENLP.NVCon.findConnection(threeCSTRDDENLP.vars.nominal.alpha);
%  
% r=threeCSTRDDENLP.vars.critical.r.values/norm(threeCSTRDDENLP.vars.critical.r.values,2);
% newNominal=threeCSTRDDENLP.vars.critical.alpha.values+r
% 
% break



threeCSTRDDENLP.initNVCons();

threeCSTRDDENLP.concatConstraints();
threeCSTRDDENLP.concatInitPoints();


%% run optimization

threeCSTRDDENLP.runOptim();
threeCSTRDDENLP.deconstructOptimum();


%% check stability
threeCSTRDDENLP.numMinEig=-0.5;

[~,threeCSTRDDENLP]=checkStabilityPoint(threeCSTRDDENLP,'nominal');
[~,threeCSTRDDENLP]=checkStabilityPoint(threeCSTRDDENLP,'critical');

% [~,threeCSTRDDENLP]=checkStabilityAtVertices(threeCSTRDDENLP);
% 
% [~,threeCSTRDDENLP]=checkStabilityAtRandom(threeCSTRDDENLP);
