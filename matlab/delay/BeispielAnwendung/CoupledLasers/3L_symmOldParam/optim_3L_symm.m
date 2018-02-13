
%% prepare matlab
close all
clear 
clc

%% define system to begin with
param=[0.2;0.2;0.2];
xGue=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'E1r'};{'E1i'};{'n1'};{'E2r'};{'E2i'};{'n2'};{'E3r'};{'E3i'};{'n3'}]);
alphaNom=VariableVector(param,9,[{'pump1'};{'pump2'};{'pump3'}]);
p=VariableVector(-0.02,12,{'omega'});

%% define system dynamics
tau=@(x,alpha,p)100/1000;

lasermod=@DDE_3L_symm;

%% collect them in an object 
aDDE_3L=DDE(@(x,xtau,alpha,p)lasermod(x,xtau,alpha,p),tau,xNom,alphaNom,p);

% clear tau
% clear alphaNom
% clear xGue
% clear param
% clear omega
% clear phi

minDist=0.01;

%% construct object ''aDDENLP''

%J=@(x)-1e7*(x(10)+x(11)+x(12)); %Maximale Pumpströme (mehr als Test gedacht) (mit p1=p2=0.1)

%  J=@(x)-1e5*(x(10)); %Maximale Pumpströme (mehr als Test gedacht) (mit p1=p2=p3=0.2)(UB=[0.5;0.4;0.4])

J=@(x)-1*abs((x(1)+1i*x(2))+(x(4)+1i*x(5))+(x(7)+1i*x(8))); %Maximale Intensität (mit p1=p2=0.1)

%J=@(x)-1e6*x(8); %Maximaler Pumpstrom (Anschmiegung an NV-Constraint) (mit p1=p2=0.2 und den UB für die Pumpströme ungleich gesetzt mit engerer Beschränkung auf dem Kanal, der nicht in der Kostenfunktion auftaucht)

%J=@(x)-1e6*(x(40)+x(71)); %Maximale abstände von der kritischen Mannigfaltigkeit (mit p1=p2=0.2)

% aDDENLP=DDENLP(J,aDDE_3L,xNom,-Inf(9,1),Inf(9,1),[0;0;0],[0.8;0.5;0.5]-sqrt(3)*minDist,-Inf,Inf);
aDDENLP=DDENLP(J,aDDE_3L,xNom,-Inf(9,1),Inf(9,1),[0;0;0],[0.8;0.5;0.5]-sqrt(3)*minDist,-Inf,Inf);
% clear aDDE_1L
% clear J
% clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-12,'TolX',1e-12);


%% initialize steady state constraints
aDDENLP.initializeStStRot();

% 
% aDDENLP.concatInitPoints();
% res = aDDENLP.stStCon.conFun(aDDENLP.initVal)
% 
% % aDDENLP.vars.nominal.x.values
% % 
% return

%% add NV Cons
%Mod-Fold Punkt (Branch2, Point 78)
xCrit1=VariableVector([0.1910;0.2940;-0.0092;0.2675;0.1776;-0.0103;0.2440;0.2395;-0.0098],0,[{'E1rcrit1'};{'E1icrit1'};{'n1crit1'};{'E2rcrit1'};{'E2icrit1'};{'n2crit1'};{'E3rcrit1'};{'E3icrit1'};{'n3crit1'}]);
alphaCrit1=VariableVector([0.1126;0.0918;0.1060],0,[{'pump1crit1'};{'pump2crit1'};{'pump3crit1'}]);
pCrit1=VariableVector(-0.0196,0,{'omegacrit1'});

aDDENLP.addNVCon('modfold',@Mod_Fold_3L_symm,@NV_Mod_Fold_3L_symm,xCrit1,alphaCrit1,pCrit1);
aDDENLP.NVCon(end).vars.w1.values=[0.502100271014960;-0.326153359036574;-2.36458521127317e-08;0.303333612451536;-0.456751674353955;-4.10886805086042e-09;0.409031372803540;-0.416628452977318;-5.15094758710284e-09];
clear alphaCrit1
clear xCrit1
clear pCrit1

%Mod-Fold Punkt (Branch2, 125)
xCrit2=VariableVector([0.1908;0.2520;-0.0103;0.0812;0.3369;-0.0093;0.1426;0.3076;-0.0097],0,[{'E1rcrit2'};{'E1icrit2'};{'n1crit2'};{'E2rcrit2'};{'E2icrit2'};{'n2crit2'};{'E3rcrit2'};{'E3icrit2'};{'n3crit2'}]);
alphaCrit2=VariableVector([0.0886;0.1097;0.1041],0,[{'pump1crit2'};{'pump2crit2'};{'pump3crit2'}]);
pCrit2=VariableVector(-0.0196,0,{'omegacrit2'});

aDDENLP.addNVCon('modfold',@Mod_Fold_3L_symm,@NV_Mod_Fold_3L_symm,xCrit2,alphaCrit2,pCrit2);
aDDENLP.NVCon(end).vars.w1.values=[0.435378614717491;-0.329732170718637;-2.99066510764544e-09;0.582100875449745;-0.140250593567891;-1.95302022478507e-08;0.531510442638750;-0.246388207125076;-4.28420307736682e-09];
clear alphacrit2
clear xcrit2
clear pcrit2

%Mod-Fold Punkt (Branch3, 137)
xCrit3=VariableVector([0.1910;0.2389;-0.0101;0.1910;0.2389;-0.0101;0.1102;0.3036;-0.0093],0,[{'E1rcrit3'};{'E1icrit3'};{'n1crit3'};{'E2rcrit3'};{'E2icrit3'};{'n2crit3'};{'E3rcrit3'};{'E3icrit3'};{'n3crit3'}]);
alphaCrit3=VariableVector([0.0826;0.0826;0.0941],0,[{'pump1crit3'};{'pump2crit3'};{'pump3crit3'}]);
pCrit3=VariableVector(-0.0197,0,{'omegacrit3'});

aDDENLP.addNVCon('modfold',@Mod_Fold_3L_symm,@NV_Mod_Fold_3L_symm,xCrit3,alphaCrit3,pCrit3);
aDDENLP.NVCon(end).vars.w1.values=[0.442605142484185;-0.353800545572045;-5.93888307558178e-09;0.442605142484183;-0.353800545572043;-5.93888306283757e-09;0.562286468057200;-0.204170642940452;-2.77089678907535e-08];
clear alphacrit2
clear xcrit2
clear pcrit2


%% initatialize constraints and prepare optimization

% aDDENLP.NVCon(1).findManifoldPoint(aDDENLP.NVCon(1).vars);
% aDDENLP.NVCon(1).findClosestCriticalPoint(alphaNom);
% aDDENLP.NVCon(1).findNormalVector(alphaNom);
% aDDENLP.NVCon(1).findConnection(alphaNom);


aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();



% openfig('Pump_3L_krit_Man.fig')

% figure(1)
% hold on
% plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'ks')
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs')
% plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs')
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),0.01*aDDENLP.vars.critical(1,1).r.values(1),0.01*aDDENLP.vars.critical(1,1).r.values(2),0.01*aDDENLP.vars.critical(1,1).r.values(3))
% quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),0.01*aDDENLP.vars.critical(2,1).r.values(1),0.01*aDDENLP.vars.critical(2,1).r.values(2),0.01*aDDENLP.vars.critical(2,1).r.values(3))
% quiver3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),0.01*aDDENLP.vars.critical(3,1).r.values(1),0.01*aDDENLP.vars.critical(3,1).r.values(2),0.01*aDDENLP.vars.critical(3,1).r.values(3))
% 

%% run optimization
aDDENLP.maxAllowedRealPart=-1;
aDDENLP.minDist=minDist;

% res = aDDENLP.stStCon.conFun(aDDENLP.initVal)

% ceq = aDDENLP.allNLEqConstraints(aDDENLP.initVal)

aDDENLP.runOptim;

aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-15,'TolX',1e-12);
aDDENLP.deconstructOptimum();
aDDENLP.concatInitPoints();
aDDENLP.runOptim;

aDDENLP.deconstructOptimum();

figure(1)
hold on
plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.alpha.values(3),'ks')

figure(1)
hold on
box on
grid on
axis equal
plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).alpha.values(3),'rs')
plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).alpha.values(3),'rs')
plot3(aDDENLP.vars.critical(3,1).alpha.values(1),aDDENLP.vars.critical(3,1).alpha.values(2),aDDENLP.vars.critical(3,1).alpha.values(3),'rs')

doubleMani1 = cross(aDDENLP.vars.critical(2,1).r.values,aDDENLP.vars.critical(3,1).r.values);
doubleMani2 = -cross(aDDENLP.vars.critical(1,1).r.values,aDDENLP.vars.critical(3,1).r.values);
doubleMani3 = cross(aDDENLP.vars.critical(1,1).r.values,aDDENLP.vars.critical(2,1).r.values);

doubleMani1 = doubleMani1/norm(doubleMani1,2);
doubleMani2 = doubleMani2/norm(doubleMani2,2);
doubleMani3 = doubleMani3/norm(doubleMani3,2);

mani1 = [[0;0;0],doubleMani2,doubleMani3];
mani2 = [[0;0;0],doubleMani1,doubleMani3];
mani3 = [[0;0;0],doubleMani1,doubleMani2];

patch(mani1(1,:),mani1(2,:),mani1(3,:),'b','FaceAlpha',.5,'EdgeColor','none')
patch(mani2(1,:),mani2(2,:),mani2(3,:),'b','FaceAlpha',.5,'EdgeColor','none')
patch(mani3(1,:),mani3(2,:),mani3(3,:),'b','FaceAlpha',.5,'EdgeColor','none')

plot3([0, doubleMani1(1)],[0, doubleMani1(2)],[0, doubleMani1(3)],'k')
plot3([0, doubleMani2(1)],[0, doubleMani2(2)],[0, doubleMani2(3)],'k')
plot3([0, doubleMani3(1)],[0, doubleMani3(2)],[0, doubleMani3(3)],'k')


xlim([0.4 0.8])
ylim([0.4 0.8])
zlim([0.4 0.8])

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

omegaDisturbance=1000*linspace(-0.5*omega,0.5*omega,3)';
absDisturbance=linspace(1,1.2,3)';

phi=0.55*2*pi;

c=exp(1i*phi);

for ii = 1:3
    E(ii)=c*(point.x.values(ii*3-2)+1i*point.x.values(3*ii-1));
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


 options=ddeset('MaxStep',0.01);
 sol=ddesd(aDDENLP,point,history,[0 2],options);
 
 figure(2);clf;
 plot(sol.x,sol.y(1:3:end,:)); 
matlab2tikz('3LsymRotSim.tex')
 
 figure(3);clf;
 A=sol.y(1:3:end,:)+1i*sol.y(2:3:end,:);
 E=A.*repmat(exp(1i*omega*sol.x*1000),3,1);
 plot(sol.x,real(E))
matlab2tikz('3LsymFixSim.tex')

