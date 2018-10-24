
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

% Startwert aus der Optimierung mit anderen Gleichungen
% param=[0.594273236532978;0.606081751647561;0.382679491924311];
% xGue=[0.707662493139026;0.331045015674947;-0.010000000000000;0.699834091868197;0.363823106221385;-0.009896375472628;0.631177784149722;0.007922629166356;-0.011275835420330];

% param=[0.5;0.4;0.35];
% xGue=[0.158150595276865;0.432495353527473;-0.00999999999999801;0.252658341977993;0.303152396987253;-0.0109845056933564;0.252658341980933;0.303152396984752;-0.0109845056933847];

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

J=@(x)-30e1*abs((x(1)+1i*x(2))+(x(4)+1i*x(5))+(x(7)+1i*x(8)))^2; %Maximale Intensität (mit p1=p2=0.1)

aDDENLP=DDENLP(J,aDDE_3L,xNom,...
    [-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf;-Inf],...
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


% Mod-Hopf
xCrit1=VariableVector([0.179745208852687;0.296078743211786;-0.00500000000000000;0.179745208852688;0.296078743211785;-0.00500000000000001;0.230481206777599;0.152647658110947;-0.00566738673785834],0,[{'E1rcrit3'};{'E1icrit3'};{'n1crit3'};{'E2rcrit3'};{'E2icrit3'};{'n2crit3'};{'E3rcrit3'};{'E3icrit3'};{'n3crit3'}]);
alphaCrit1=VariableVector([0.113771252664493;0.113771252664493;0.0698892712723909],0,[{'pump1crit3'};{'pump2crit3'};{'pump3crit3'}]);
pCrit1=VariableVector(-0.0200000000000000,0,{'omegacrit3'});

aDDENLP.addNVCon('modhopf',@Mod_Hopf_3L_SETUP_1,@NV_Mod_Hopf_3L_SETUP_1,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.omega.values=28.802287251711665;
aDDENLP.NVCon(end).vars.w1.values=[-0.000500125722836925- 1.01864366297517e-13i;0.000504383652331200 + 1.02781354257418e-13i;-1.53683297937919e-06 + 1.61939958350911e-15i;2.77822399212342e-05 + 1.52681084191203e-13i;-2.76591072242404e-05 - 1.54167788695597e-13i;7.58183255391134e-06 + 4.74216609600320e-15i;-0.0851091237755987 - 4.60413473787358e-14i;-0.0288044059060574 + 1.56669538457187e-13i;0.0249533489218718 + 2.09701251318971e-16i];
aDDENLP.NVCon(end).vars.w2.values=[-0.000206247903144366 + 2.22920241066234e-14i;0.000206448261575178 - 2.27215129616886e-14i;1.69443234398295e-05 + 2.85559071153774e-15i;0.000185830355098729 + 1.66808610250407e-13i;-0.000187795330143035 - 1.68084254016739e-13i;2.70356615987553e-07 - 5.01352556309493e-15i;0.318108137531340 - 3.91747048632784e-14i;-0.943447571035896 + 3.35315257696474e-14i;-0.00407799152273783 + 4.27224547493030e-15i];
clear alphaCrit1
clear xCrit1
clear pCrit1


% Mod-Hopf
xCrit1=VariableVector([0.0481993460120035;0.405218965350386;-0.00500000000000000;0.172609612890016;0.276219407092834;-0.00566724517419908;0.0475581323195376;0.405461593871524;-0.00499796348065467],0,[{'E1rcrit4'};{'E1icrit4'};{'n1crit4'};{'E2rcrit4'};{'E2icrit4'};{'n2crit4'};{'E3rcrit4'};{'E3icrit4'};{'n3crit4'}]);
alphaCrit1=VariableVector([0.159860330967266;0.0992215040144732;0.159996986589572],0,[{'pump1crit4'};{'pump2crit4'};{'pump3crit4'}]);
pCrit1=VariableVector(-0.0200000000000000,0,{'omegacrit4'});

aDDENLP.addNVCon('modhopf',@Mod_Hopf_3L_SETUP_1,@NV_Mod_Hopf_3L_SETUP_1,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.omega.values=33.548686577848656;
aDDENLP.NVCon(end).vars.w1.values=[0.000211370539298167 - 2.61625857212135e-11i;-8.03645115827106e-05 + 9.96496712356446e-12i;5.21954908663972e-06 - 3.82016772629373e-13i;0.0564922410314140 + 4.45070487056059e-10i;0.0534207237171074 - 4.43502997481319e-10i;-0.0247453741012018 - 7.63772716626017e-12i;5.91733477053411e-06 + 7.41743451262085e-12i;-2.35734327081444e-06 - 2.82223903082954e-12i;-3.35835613493908e-06 + 1.60473448499753e-12i];
aDDENLP.NVCon(end).vars.w2.values=[0.000315140908886033 - 2.61605914884881e-11i;-0.000118742842840138 + 9.89027100463734e-12i;-6.23274587415935e-06 + 6.95009739741382e-13i;-0.683982215524306 - 1.20137712314790e-10i;0.724912608484712 + 1.95779623813987e-10i;0.00346505392614992 - 1.47353813425001e-11i;-0.000111874689768821 + 7.10793352925624e-11i;4.22721643390300e-05 - 2.68047688436548e-11i;-7.66274202976462e-07 - 2.68164190796746e-13i];
clear alphaCrit1
clear xCrit1
clear pCrit1


%% initatialize constraints and prepare optimization

aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();

%% plot

openfig('crit_Manifold_3L_asymm_Old_Param_large_clear.fig')

figure(1)
hold on
plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'ks','LineWidth',2)

figure(1)
hold on
plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs','LineWidth',2)
plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs','LineWidth',2)
plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs','LineWidth',2)
plot3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),'rs','LineWidth',2)


figure(1)
hold on
quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.055*aDDENLP.vars.critical(1,1).r.values(1),0.055*aDDENLP.vars.critical(1,1).r.values(2),0.055*aDDENLP.vars.critical(1,1).r.values(3),'LineWidth',2)
quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.055*aDDENLP.vars.critical(2,1).r.values(1),0.055*aDDENLP.vars.critical(2,1).r.values(2),0.055*aDDENLP.vars.critical(2,1).r.values(3),'LineWidth',2)
quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),0.38*aDDENLP.vars.critical(3,1).r.values(1),0.38*aDDENLP.vars.critical(3,1).r.values(2),0.38*aDDENLP.vars.critical(3,1).r.values(3),'LineWidth',2)
quiver3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),0.38*aDDENLP.vars.critical(4,1).r.values(1),0.38*aDDENLP.vars.critical(4,1).r.values(2),0.38*aDDENLP.vars.critical(4,1).r.values(3),'LineWidth',2)

xlim([0 0.8])
ylim([0 0.8])
zlim([0 0.8])


%% run optimization
aDDENLP.maxAllowedRealPart=-1;
aDDENLP.minDist=minDist;

aDDENLP.runOptim;

aDDENLP.deconstructOptimum();

aDDENLP.aCostFunction=@(x)-abs((x(1)+1i*x(2))+(x(4)+1i*x(5))+(x(7)+1i*x(8)))^2;
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-12,'TolX',1e-12);
aDDENLP.deconstructOptimum();
aDDENLP.concatInitPoints();
aDDENLP.runOptim;

aDDENLP.deconstructOptimum();

%% plot
figure(1)
hold on
plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'ks','LineWidth',2)

figure(1)
hold on
plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs','LineWidth',2)
plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs','LineWidth',2)
plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs','LineWidth',2)
plot3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),'rs','LineWidth',2)
figure(1)
hold on
quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.02*aDDENLP.vars.critical(1,1).r.values(1),0.02*aDDENLP.vars.critical(1,1).r.values(2),0.02*aDDENLP.vars.critical(1,1).r.values(3),'LineWidth',2)
quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.02*aDDENLP.vars.critical(2,1).r.values(1),0.02*aDDENLP.vars.critical(2,1).r.values(2),0.02*aDDENLP.vars.critical(2,1).r.values(3),'LineWidth',2)
quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),7.7*aDDENLP.vars.critical(3,1).r.values(1),7.7*aDDENLP.vars.critical(3,1).r.values(2),7.7*aDDENLP.vars.critical(3,1).r.values(3),'LineWidth',2)
quiver3(aDDENLP.vars.critical(4,1).alpha.values(1),aDDENLP.vars.critical(4,1).alpha.values(2),aDDENLP.vars.critical(4,1).alpha.values(3),7.7*aDDENLP.vars.critical(4,1).r.values(1),7.7*aDDENLP.vars.critical(4,1).r.values(2),7.7*aDDENLP.vars.critical(4,1).r.values(3),'LineWidth',2)
box on

xlim([0 0.8])
ylim([0 0.8])
zlim([0 0.8])

%% 

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

%% manually put the nominal point to the right mode

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

 options=ddeset('RelTol',1e-12,'AbsTol',1e-9);
 sol=ddesd(aDDENLP,point,history,[0 2],options);
 
 aDDENLP.vars.nominal.x.values=sol.y(:,end)
% 



%% check stability
%aDDENLP.vars.nominal.x.values=aDDENLP.vars.nominal.x.values.*[1,1,1,-1,-1,1,-1,-1,1]';
[~]=aDDENLP.checkStabilityPoint('nominal');
[~]=aDDENLP.checkStabilityPoint('critical');
% 
% [~]=aDDENLP.checkStabilityAtVertices();
% 
% [~]=aDDENLP.checkStabilityAtRandom();


%% run simulation
 
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


 options=ddeset('RelTol',1e-12,'AbsTol',1e-9);
 sol=ddesd(aDDENLP,point,history,[0 2],options);
 
 
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

%% 
