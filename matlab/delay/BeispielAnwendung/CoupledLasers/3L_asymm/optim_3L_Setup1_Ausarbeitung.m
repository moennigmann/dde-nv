
%% prepare matlab
close all
clear
clc

%% define system to begin with
% param=[0.22;0.22;0.22];
% xGue=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01];

% Startwert aus findstart
param=[0.199943165310036;0.143042379378885;0.143042379378884]*2;
xGue=[0.158150595276865;0.432495353527473;-0.00999999999999801;0.252658341977993;0.303152396987253;-0.0109845056933564;0.252658341980933;0.303152396984752;-0.0109845056933847];

param=[0.5;0.4;0.35];
xGue=[0.158150595276865;0.432495353527473;-0.00999999999999801;0.252658341977993;0.303152396987253;-0.0109845056933564;0.252658341980933;0.303152396984752;-0.0109845056933847];

% Startwert aus erster Optimierung mit Zustandsgrößen aus fsolve
% (Kostenfunktion)
% param=[0.499999986535985;0.506111305451097;0.506111305451122];
% xGue=[0.427666820986429;0.576413571473567;-0.0100000000000000;0.406994048466545;0.599540183469653;-0.00989617758939340;0.406994048466545;0.599540183469653;-0.00989617758939340];

% Zustände nach Einschwingvorgang
% xGue=[-0.0710637457134606;0.714213865449574;-0.0100000000018390;-0.0892431237220643;0.716425249839996;-0.00993846030598379;-0.0892431237065467;0.716425249837994;-0.00993846030750198];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'E1r'};{'E1i'};{'n1'};{'E2r'};{'E2i'};{'n2'};{'E3r'};{'E3i'};{'n3'}]);
alphaNom=VariableVector(param,9,[{'pump1'};{'pump2'};{'pump3'}]);
p=VariableVector(-0.02,12,{'omega'});

minDist = 0.01;

%% define system dynamics
tau=@(x,alpha,p)100/1000;

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

J=@(x)-10000*abs((x(1)+1i*x(2))+(x(4)+1i*x(5))+(x(7)+1i*x(8)))^2; %Maximale Intensität (mit p1=p2=0.1)

aDDENLP=DDENLP(J,aDDE_3L,xNom,...
    [0;0;-Inf;0;0;-Inf;0;0;-Inf],...
    [Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf;Inf],...
    [0;0;0]+sqrt(3)*minDist,...
    [0.8;0.7;0.4]-sqrt(3)*minDist,...
    -Inf,Inf);

% clear aDDE_1L
% clear J
% clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-16,'TolX',1e-12);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',2000000,'display','iter','TolFun',1e-12,'TolX',1e-12,'ScaleProblem','none');
%active-set
%% initialize steady state constraints
aDDENLP.initializeStStRot();

%% add NV Cons
% Mod-Fold Punkt (Branch3, Point 125)
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
aDDENLP.NVCon(end).vars.omega.values=0.0213956378669933e3;
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
aDDENLP.NVCon(end).vars.omega.values=0.0226416709337338e3;
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


%% run optimization
aDDENLP.maxAllowedRealPart=-1;
aDDENLP.minDist=minDist;

aDDENLP.runOptim;

aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-12,'TolX',1e-12);
aDDENLP.deconstructOptimum();
aDDENLP.concatInitPoints();
aDDENLP.runOptim;

aDDENLP.deconstructOptimum();


disp('alphaNom =')
disp(aDDENLP.vars.nominal.alpha.values);
disp('l =');
disp([aDDENLP.vars.critical(1).l.values;aDDENLP.vars.critical(2).l.values;aDDENLP.vars.critical(3).l.values;aDDENLP.vars.critical(4).l.values]);

% %% plot results
% 
% figure(1);clf;
% hold on
% axis equal
% box on
% grid on
% 
% plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'x')
% plot3(aDDENLP.vars.critical(2).alpha.values(1),aDDENLP.vars.critical(2).alpha.values(2),aDDENLP.vars.critical(2).alpha.values(3),'xb')
% plot3(aDDENLP.vars.critical(4).alpha.values(1),aDDENLP.vars.critical(4).alpha.values(2),aDDENLP.vars.critical(4).alpha.values(3),'xk')
% 
% 
% 
% plotcube(2*minDist*ones(1,3),aDDENLP.vars.nominal.alpha.values'-minDist,0.2,[1 0 0])
% 
% [x,y,z]=sphere(20);
% surf(minDist*sqrt(3)*x+aDDENLP.vars.nominal.alpha.values(1),...
%     minDist*sqrt(3)*y+aDDENLP.vars.nominal.alpha.values(2),...
%     minDist*sqrt(3)*z+aDDENLP.vars.nominal.alpha.values(3),...
%     'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');
% 
% xlim([0.4,0.69])
% ylim([0.4,0.69])
% zlim([0.2,0.5])
% 
% fill3([0.8,0.8,0.8,0.8],[0,0,0.8,0.8],[0,0.8,0.8,0],[0.2,0.2,0.2])
% fill3([0,0,0.8,0.8],[0.7,0.7,0.7,0.7],[0,0.8,0.8,0],[0.2,0.2,0.2])
% fill3([0,0,0.8,0.8],[0,0.8,0.8,0],[0.4,0.4,0.4,0.4],[0.2,0.2,0.2])



%% check stability

[~]=aDDENLP.checkStabilityPoint('nominal');
% [~]=aDDENLP.checkStabilityPoint('critical');
% 
% [~]=aDDENLP.checkStabilityAtVertices();
% 
% [~]=aDDENLP.checkStabilityAtRandom();

return

%% run simulation

 point=aDDENLP.vars.nominal; 
 
%  history=point.x.values+0.001*randn(30,1);
omega=point.p.values(1);
history=@(t)[];

phi=[0,2*pi/6,-pi];

omegaDisturbance=[0;-0.5*omega;0.3*omega];
absDisturbance=[1,1.3,1.2];%linspace(1,1,3)';

for ii = 1:3
    E(ii)=exp(1i*phi(ii))*(point.x.values(ii*3-2)+1i*point.x.values(3*ii-1));
    n(ii)=point.x.values(ii*3);
    history=@(t)[history(t);real(E(ii)*absDisturbance(ii)*exp(+1i*omegaDisturbance(ii)*t));imag(E(ii)*absDisturbance(ii)*exp(+1i*omegaDisturbance(ii)*t));n(ii)+t*0];
end


% for ii=1:10
%    Ai=history(3*(ii-1)+1)+1i*history(3*(ii-1)+2)*exp(1i*0.01*(rand(1)-0.5)*pi);
%    history(3*(ii-1)+1)=real(Ai);
%    history(3*(ii-1)+2)=imag(Ai);
% end
% 
% history=repmat([0;0;(point.x.values(3)+point.x.values(18))/2],10,1)+0.001*randn(30,1);


 options=ddeset('MaxStep',1,'RelTol',1e-12,'AbsTol',1e-9);
 sol=ddesd(aDDENLP,point,history,[0 2000],options);
 
 figure(2);clf;
 plot(sol.x/1000,sol.y(1:3:end,:)); 
 legend('laser 1','laser 2','laser 3')
 matlab2tikz('3LasymStabRotSim.tex')
 
 
 figure(3);clf;
 A=sol.y(1:3:end,:)+1i*sol.y(2:3:end,:);
 E=A.*repmat(exp(1i*omega*sol.x),3,1);
 plot(sol.x/1000,real(E))
 legend('laser 1','laser 2','laser 3')
 matlab2tikz('3LasymStabFixSim.tex')

