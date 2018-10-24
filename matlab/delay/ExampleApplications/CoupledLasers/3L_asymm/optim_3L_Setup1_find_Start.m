
%% prepare matlab
close all
clear all
clc

%% define system to begin with
param=[0.2;0.2;0.2];
xGue=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'E1r'};{'E1i'};{'n1'};{'E2r'};{'E2i'};{'n2'};{'E3r'};{'E3i'};{'n3'}]);
alphaNom=VariableVector(param,9,[{'pump1'};{'pump2'};{'pump3'}]);
p=VariableVector(-0.02,12,{'omega'});

%% define system dynamics
tau=@(x,alpha,p)100;

lasermod=@DDE_3L_SETUP_1;

%% collect them in an object 
aDDE_3L=DDE(@(x,xtau,alpha,p)lasermod(x,xtau,alpha,p),tau,xNom,alphaNom,p);

% clear tau
% clear alphaNom
% clear xGue
% clear param
% clear omega
% clear phi

%% construct object ''aDDENLP''

%J=@(x)-1e3*(x(10)); %Maximale Pumpströme (mehr als Test gedacht) (mit p1=p2=0.1)

%J=@(x)-1e6*(x(10)); %Maximale Pumpströme (mehr als Test gedacht) (mit p1=p2=p3=0.2)(UB=[0.5;0.4;0.4])

%J=@(x)-1e3*abs((x(4)+1i*x(5))+(x(7)+1i*x(8))); %Maximale Intensität (mit p1=p2=0.1)

J=@(x)-1e3*(5*x(58)+5*x(103)+x(168)+x(233)); %Maximaler abstand zu constraints

%J=@(x)-1e6*(x(40)+x(71)); %Maximale abstände von der kritischen Mannigfaltigkeit (mit p1=p2=0.2)

aDDENLP=DDENLP(J,aDDE_3L,xNom,[-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf],[Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf],[0;0;0],[0.2;0.5;0.5],[-0.021],[-0.019]);
% clear aDDE_1L
% clear J
% clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-6,'TolX',1e-6);

%active-set
%% initialize steady state constraints
aDDENLP.initializeStSt();


%% add NV Cons
%Mod-Fold Punkt (Branch3, Point 125)
xCrit1=VariableVector([0.191000000000005;0.251111759568000;-0.0100000000000000;0.191000000000005;0.251111759568000;-0.0100000000000000;0.141034526944710;0.292687948139146;-0.00951476512316607],0,[{'E1rcrit1'};{'E1icrit1'};{'n1crit1'};{'E2rcrit1'};{'E2icrit1'};{'n2crit1'};{'E3rcrit1'};{'E3icrit1'};{'n3crit1'}]);
alphaCrit1=VariableVector([0.0885427346354054;0.0885427346354055;0.0950378578501752],0,[{'pump1crit1'};{'pump2crit1'};{'pump3crit1'}]);
pCrit1=VariableVector(-0.0200000000000000,0,{'omegacrit1'});

aDDENLP.addNVCon('modfold',@Mod_Fold_3L_SETUP_1,@NV_Mod_Fold_3L_SETUP_1,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.w1.values=[-0.000117954796360985;-7.44048497514761e-05;-3.24883500759214e-08;-0.000611842645172190;0.000457391263095490;4.41471711249944e-06;-0.904315586328776;0.426824499491078;0.00579343850667159];
clear alphaCrit1
clear xCrit1
clear pCrit1

%Mod-Fold Punkt (Branch2, Point 114)
xCrit1=VariableVector([0.191000000000010;0.261899880887572;-0.0100000000000000;0.138768370102688;0.303596595781310;-0.00951440354043023;0.191000000000010;0.261899880887572;-0.0100000000000000],0,[{'E1rcrit2'};{'E1icrit2'};{'n1crit2'};{'E2rcrit2'};{'E2icrit2'};{'n2crit2'};{'E3rcrit2'};{'E3icrit2'};{'n3crit2'}]);
alphaCrit1=VariableVector([0.0940218221328392;0.100852983260900;0.0940218221328393],0,[{'pump1crit2'};{'pump2crit2'};{'pump3crit2'}]);
pCrit1=VariableVector(-0.0200000000000000,0,{'omegacrit2'});

aDDENLP.addNVCon('modfold',@Mod_Fold_3L_SETUP_1,@NV_Mod_Fold_3L_SETUP_1,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.w1.values=[7.29473648197397e-05;3.27028722192795e-05;4.71781265486850e-08;-0.912708223566922;0.408572769756642;0.00563941934138202;0.000345630258742933;-0.000248072220509548;-2.28677330833367e-06];
clear alphaCrit1
clear xCrit1
clear pCrit1

%Mod-Hopf Punkt (Branch2, Point 86)
xCrit1=VariableVector([0.190999999999984;0.300389769642497;-0.0100000000000000;0.239095829192106;0.144995370604579;-0.0114109079903168;0.190999999999984;0.300389769642497;-0.0100000000000000],0,[{'E1rcrit3'};{'E1icrit3'};{'n1crit3'};{'E2rcrit3'};{'E2icrit3'};{'n2crit3'};{'E3rcrit3'};{'E3icrit3'};{'n3crit3'}]);
alphaCrit1=VariableVector([0.115447863568808;0.0658873407499950;0.115447863568808],0,[{'pump1crit3'};{'pump2crit3'};{'pump3crit3'}]);
pCrit1=VariableVector(-0.0200000000000000,0,{'omegacrit3'});

aDDENLP.addNVCon('modhopf',@Mod_Hopf_3L_SETUP_1,@NV_Mod_Hopf_3L_SETUP_1,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.omega.values=0.0213956378669933;
aDDENLP.NVCon(end).vars.w1.values=[-0.000230206260101463 + 4.65240431575402e-15i;0.000241082143039494 - 4.91882666665753e-15i;1.10687798312290e-05 + 1.81975368316743e-16i;0.119998673283529 - 9.34677877135260e-16i;0.0337931522999705 - 2.01965265392279e-16i;-0.0361985288850473 + 9.86115091444077e-18i;7.05700680574231e-05 + 8.95037702699407e-15i;-7.39428817757922e-05 - 9.37460921611299e-15i;6.65444476584798e-06 - 4.20477787709720e-16i];
aDDENLP.NVCon(end).vars.w2.values=[0.000152998833607872 + 4.78448265711729e-15i;-0.000161778632058078 - 5.02123865012095e-15i;9.72087566066553e-06 - 3.06191670125890e-16i;-0.266498801752518 - 4.89265869719861e-16i;0.955018559613074 - 8.71491330189171e-18i;0.00810879499551736 - 7.92561040389832e-17i;0.000105028449735736 - 8.65617104004639e-15i;-0.000110606754458253 + 9.11156988969458e-15i;-2.59963687879822e-06 - 3.80654289736055e-16i];
clear alphaCrit1
clear xCrit1
clear pCrit1

%Mod-Hopf Punkt (Branch2, Point 88)
xCrit1=VariableVector([0.191000000000013;0.326363903211633;-0.0100000000000000;0.191000000000013;0.326363903211633;-0.0100000000000000;0.248459536467490;0.166036073417824;-0.0113822863473270],0,[{'E1rcrit4'};{'E1icrit4'};{'n1crit4'};{'E2rcrit4'};{'E2icrit4'};{'n2crit4'};{'E3rcrit4'};{'E3icrit4'};{'n3crit4'}]);
alphaCrit1=VariableVector([0.131564453346341;0.131564453346342;0.0769013930657235],0,[{'pump1crit4'};{'pump2crit4'};{'pump3crit4'}]);
pCrit1=VariableVector(-0.0200000000000000,0,{'omegacrit4'});

aDDENLP.addNVCon('modhopf',@Mod_Hopf_3L_SETUP_1,@NV_Mod_Hopf_3L_SETUP_1,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.omega.values=0.0226416709337338;
aDDENLP.NVCon(end).vars.w1.values=[0.000220997495238136 - 5.41476002464851e-20i;-0.000215088997583262 + 4.95654635851059e-20i;-7.67915525812985e-06 + 4.48654643589372e-20i;-5.10480105309870e-05 + 1.77993219485698e-18i;4.96269916047288e-05 - 1.73568809483612e-18i;-5.52832461453007e-06 - 7.22053197845304e-21i;-0.110189446077025 - 5.62020896310451e-20i;-0.0367557827574800 - 3.69999080383496e-21i;0.0359638679763530 + 7.83268613949092e-21i];
aDDENLP.NVCon(end).vars.w2.values=[-9.71296390333951e-05 + 9.14925707762996e-19i;9.57131982807702e-05 - 8.95119369938142e-19i;-9.32358518663881e-06 - 9.81574016118904e-21i;-9.23706654649083e-05 - 4.21836183163345e-20i;9.03659390567647e-05 + 4.52949549410834e-20i;1.66738478040298e-06 - 7.73039319425692e-20i;0.311846768410686 - 1.18572666218663e-19i;-0.942288717383010 - 3.19254461045111e-20i;-0.00757159739032051 - 6.05344583846413e-21i];
clear alphaCrit1
clear xCrit1
clear pCrit1

%% initatialize constraints and prepare optimization

% aDDENLP.NVCon(1).findManifoldPoint(aDDENLP.NVCon(1).vars);
% aDDENLP.NVCon(1).findClosestCriticalPoint(alphaNom);
% aDDENLP.NVCon(1).findNormalVector(alphaNom);
% aDDENLP.NVCon(1).findConnection(alphaNom);

aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();

openfig('Guete_p1_p2.fig')
openfig('Guete_p2_p3.fig')

figure(1)
hold on
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'ks')

figure(2)
hold on
plot(aDDENLP.vars.nominal.alpha.values(3),aDDENLP.vars.nominal.alpha.values(2),'ks')


% figure(1)
% hold on
% plot(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),'rs')
% plot(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),'rs')
% 
% figure(2)
% hold on
% plot(aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs')
% plot(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(3),'rs')
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.01*aDDENLP.vars.critical(1,1).r.values(1),0.01*aDDENLP.vars.critical(1,1).r.values(2),0.01*aDDENLP.vars.critical(1,1).r.values(3))
% quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.01*aDDENLP.vars.critical(2,1).r.values(1),0.01*aDDENLP.vars.critical(2,1).r.values(2),0.01*aDDENLP.vars.critical(2,1).r.values(3))
% quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),0.5*aDDENLP.vars.critical(3,1).r.values(1),0.5*aDDENLP.vars.critical(3,1).r.values(2),0.5*aDDENLP.vars.critical(3,1).r.values(3))
% quiver3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),0.5*aDDENLP.vars.critical(4,1).r.values(1),0.5*aDDENLP.vars.critical(4,1).r.values(2),0.5*aDDENLP.vars.critical(4,1).r.values(3))
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),'rs')
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.01*aDDENLP.vars.critical(1,1).r.values(1),0.01*aDDENLP.vars.critical(1,1).r.values(2),0.01*aDDENLP.vars.critical(1,1).r.values(3))
% quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.01*aDDENLP.vars.critical(2,1).r.values(1),0.01*aDDENLP.vars.critical(2,1).r.values(2),0.01*aDDENLP.vars.critical(2,1).r.values(3))
% quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),0.5*aDDENLP.vars.critical(3,1).r.values(1),0.5*aDDENLP.vars.critical(3,1).r.values(2),0.5*aDDENLP.vars.critical(3,1).r.values(3))
% quiver3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),0.5*aDDENLP.vars.critical(4,1).r.values(1),0.5*aDDENLP.vars.critical(4,1).r.values(2),0.5*aDDENLP.vars.critical(4,1).r.values(3))
% 


%% run optimization
aDDENLP.maxAllowedRealPart=-0.001;
aDDENLP.minDist=0.01;

aDDENLP.runOptim;
hier = aDDENLP.deconstructOptimum();

figure(1)
hold on
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'ks')

figure(2)
hold on
plot(aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'ks')


% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'ks')
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),'rs')
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.1*aDDENLP.vars.critical(1,1).r.values(1),0.1*aDDENLP.vars.critical(1,1).r.values(2),0.1*aDDENLP.vars.critical(1,1).r.values(3))
% quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.1*aDDENLP.vars.critical(2,1).r.values(1),0.1*aDDENLP.vars.critical(2,1).r.values(2),0.1*aDDENLP.vars.critical(2,1).r.values(3))
% quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),3*aDDENLP.vars.critical(3,1).r.values(1),3*aDDENLP.vars.critical(3,1).r.values(2),3*aDDENLP.vars.critical(3,1).r.values(3))
% quiver3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),3*aDDENLP.vars.critical(4,1).r.values(1),3*aDDENLP.vars.critical(4,1).r.values(2),3*aDDENLP.vars.critical(4,1).r.values(3))
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),'rs')
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.1*aDDENLP.vars.critical(1,1).r.values(1),0.1*aDDENLP.vars.critical(1,1).r.values(2),0.1*aDDENLP.vars.critical(1,1).r.values(3))
% quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.1*aDDENLP.vars.critical(2,1).r.values(1),0.1*aDDENLP.vars.critical(2,1).r.values(2),0.1*aDDENLP.vars.critical(2,1).r.values(3))
% quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),3*aDDENLP.vars.critical(3,1).r.values(1),3*aDDENLP.vars.critical(3,1).r.values(2),3*aDDENLP.vars.critical(3,1).r.values(3))
% quiver3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),3*aDDENLP.vars.critical(4,1).r.values(1),3*aDDENLP.vars.critical(4,1).r.values(2),3*aDDENLP.vars.critical(4,1).r.values(3))
% 
% 
% 

%% run simulation

%  point=aDDENLP.vars.nominal;
%  
%  history=abs(point.x.values);
%  options=ddeset();
%  sol=ddesd(aDDENLP,point,history,[0 5000],options);
%  
%  figure(2);clf;
%  plot(sol.x,sol.y);

%% check stability

% [~]=aDDENLP.checkStabilityPoint('nominal');
% [~]=aDDENLP.checkStabilityPoint('critical');
% 
% [~]=aDDENLP.checkStabilityAtVertices();
% 
% [~]=aDDENLP.checkStabilityAtRandom();
