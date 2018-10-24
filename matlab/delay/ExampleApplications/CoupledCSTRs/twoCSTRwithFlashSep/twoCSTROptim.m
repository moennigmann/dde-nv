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
Fr=5; % m^3/h
beta=15*10;

T01=400; %K
F01=5;% m^3/h
Tc1=480; %K
V1=1*1000; %m^3

T02=405; %K
F02=5;% m^3/h
Tc2=485; %K
V2=2*1000; %m^3

TcF=480; %K
VF=1*1000; %m^3

param=[Fr, beta, T01, F01, Tc1, V1, T02, F02, Tc2, V2, TcF, VF]';
xGue=[0.070465268853638;0.476596695860939;5.074024028147752e+02;0.069117014670805;0.457104951396088;5.052625397149701e+02;0.032816753806440;0.424051142574560;4.800874136322317e+02];
%% name states and parameters
xNom=VariableVector(xGue,0,[{'cA1'};{'cB1'};{'T1'};{'cA2'};{'cB2'};{'T2'};{'cAF'};{'cBF'};{'TF'}]);
alphaNom=VariableVector(param,9,[{'Fr'}; {'beta'};{'T01'};{'F01'};{'Tc1'};{'V1'};{'T02'};{'F02'};{'Tc2'};{'V2'};{'TcF'};{'VF'}]);
p=VariableVector([],9+12,[]);



clear Fr T01 F01 Tc1 V1 T02 F02 Tc2 V2 TcF VF beta xGue param


%% define system dynamics
tau=@(x,alpha,p)[5*350/(alpha(1));
    5*350/(alpha(1)+alpha(4));
    5*350/(alpha(1)+alpha(4)+alpha(8))];

cstrModel=@(x,xtau,alpha,p)twoCSTRoneFlashSepDDE(x,xtau,alpha,p);



%% collect them in an object 
cstrDDE=DDE(@(x,xtau,alpha,p)0.01*cstrModel(x,xtau,alpha,p),tau,xNom,alphaNom,p);


clear cstrModel
clear tau
clear alphaNom


%% construct object ''aDDENLP''

J=@(x)100*returnOnInvest2(x);

setNCSTRBoundaries


twoCSTRDDENLP=DDENLP(J,cstrDDE,xNom,...
    [Cmin;Cmin;Tmin;Cmin;Cmin;Tmin;Cmin;Cmin;Tmin],...
    [Cmax;Cmax;Tmax;Cmax;Cmax;Tmax;Cmax;Cmax;Tmax],...
    [Frmin, betamin, Tmin, F0min, Tcmin, Vmin, Tmin, F0min, Tcmin, Vmin, Tcmin, Vmin]'+minDist*sqrt(12),...
    [Frmax, betamax, T0max, F0max, Tmax, Vmax, T0max, F0max, Tmax, Vmax, Tmax, Vmax]'-minDist*sqrt(12),...
    -Inf(0,1),...
    Inf(0,1));

clear cstrDDE J xNom Cmin Tmin Cmax Tmax Frmin betamin F0min Vmin Frmax betamax Tmax F0max Vmax



twoCSTRDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',100000,'MaxFunEvals',2000000,'display','off','TolFun',1e-7,'TolX',1e-7);
twoCSTRDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-9);
twoCSTRDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-9,'TolX',1e-9,'TolCon',1e-6);


%% initialize steady state constraints
twoCSTRDDENLP.initializeStSt();


%% add NV Cons

% first
alphaCrit1=VariableVector([5 119.809004809229 400 8.29919026944543 480 1000 400 5 480 1000 480 1000]',0,{'alphacrit1'});
xCrit1=VariableVector([0.079170911838233;0.474090456939760;5.231186142094912e+02;0.069535797060128;0.460759938311193;5.060187356472061e+02;0.036981097893851;0.436182179170303;4.801099227099807e+02],0,{'xcrit1'});

twoCSTRDDENLP.addNVCon('hopf',@twoCSTRoneFlashSepHopfMani,@twoCSTRoneFlashSepHopfManiNV,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1

twoCSTRDDENLP.NVCon.vars.omega.values=0.1343721963052;
twoCSTRDDENLP.NVCon.vars.w1.values=real([-7.819987534297238e-05 + 3.460604469251216e-04i;7.618235205292739e-05 - 1.074230213547849e-04i;0.999842086609801 - 0.010575510422037i;8.406094691683397e-06 + 6.909989412539884e-06i;-1.803731083489705e-06 - 4.464266896916172e-06i;0.014266765095470 - 0.000522384456407i;2.789837867784241e-07 - 3.013531185467610e-07i;-1.764504724490689e-07 + 5.740329800094467e-08i;5.922694062575986e-05 - 9.226400684842233e-06i]);
twoCSTRDDENLP.NVCon.vars.w2.values=imag([-7.819987534297238e-05 + 3.460604469251216e-04i;7.618235205292739e-05 - 1.074230213547849e-04i;0.999842086609801 - 0.010575510422037i;8.406094691683397e-06 + 6.909989412539884e-06i;-1.803731083489705e-06 - 4.464266896916172e-06i;0.014266765095470 - 0.000522384456407i;2.789837867784241e-07 - 3.013531185467610e-07i;-1.764504724490689e-07 + 5.740329800094467e-08i;5.922694062575986e-05 - 9.226400684842233e-06i]);

% second 
alphaCrit2=VariableVector([22.1601213359225 168.682361101347 400.021282729215 2.27841619729643 490.907117874265 1001.08552704709 399.929237964557 29.8743994643126 451.162608418820 1001.85565115951 480.017815727392 1001.67368297621]',0,{'alphacrit2'});
xCrit2=VariableVector([0.083708826644944;0.541731318027979;5.358063493806366e+02;0.312310746953949;0.485779053884835;4.910082755357364e+02;0.187872433111242;0.541071115653070;4.801151352648811e+02],0,{'xcrit2'});

twoCSTRDDENLP.addNVCon('hopf',@twoCSTRoneFlashSepHopfMani,@twoCSTRoneFlashSepHopfManiNV,xCrit2,alphaCrit2,p);
clear alphaCrit2
clear xCrit2

twoCSTRDDENLP.NVCon(2).vars.omega.values=0.1343721963052;
twoCSTRDDENLP.NVCon(2).vars.w1.values=real([-4.537804476178019e-06 + 6.056157450158904e-06i;7.031103721776286e-07 + 3.157524524377460e-06i;0.028336452684666 + 0.009799324327360i;-1.613734475939913e-04 + 7.510348494881957e-04i;1.537206892741151e-04 - 6.362720553358742e-04i;0.999510940681738 + 0.000000039205764i;8.169443439284826e-05 + 2.980704105737675e-05i;-6.904975211786184e-05 - 2.719158298653648e-05i;0.008793220730219 - 0.000738726302095i]);
twoCSTRDDENLP.NVCon(2).vars.w2.values=imag([-4.537804476178019e-06 + 6.056157450158904e-06i;7.031103721776286e-07 + 3.157524524377460e-06i;0.028336452684666 + 0.009799324327360i;-1.613734475939913e-04 + 7.510348494881957e-04i;1.537206892741151e-04 - 6.362720553358742e-04i;0.999510940681738 + 0.000000039205764i;8.169443439284826e-05 + 2.980704105737675e-05i;-6.904975211786184e-05 - 2.719158298653648e-05i;0.008793220730219 - 0.000738726302095i]);

%% initatialize constraints and prepare optimization

twoCSTRDDENLP.initNVCons();

%alphaCritClosest1=
% 4.8642865      149.16322      400.00166      12.037311      480.73989      999.95926      404.99987       5.567524      484.80748      1999.9975      480.00071      999.99999
%alphaCritClosest2=
% 6.0102816      148.32749      399.99862      7.1791801      479.14273      999.98204      405.00135      15.796538      486.37059      1999.9526      479.99868      1000.0003


twoCSTRDDENLP.concatConstraints();
twoCSTRDDENLP.concatInitPoints();


%% run optimization

twoCSTRDDENLP.runOptim();
twoCSTRDDENLP.deconstructOptimum();


%% check stability
twoCSTRDDENLP.numMinEig=-0.5;

maxEigN=twoCSTRDDENLP.checkStabilityPoint('nominal');
maxEigC=twoCSTRDDENLP.checkStabilityPoint('critical');
% 
% [~,oneCSTRDDENLP]=checkStabilityAtVertices(oneCSTRDDENLP);
% [~,oneCSTRDDENLP]=checkStabilityAtRandom(oneCSTRDDENLP);
