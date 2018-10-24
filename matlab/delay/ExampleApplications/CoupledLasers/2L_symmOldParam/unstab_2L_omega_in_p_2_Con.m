
%% prepare matlab
close all
clear 
clc

%% define system to begin with

param = [0.45;0.45];
xGue = [0.0407031314951728;0.668830467436600;0.000696388699168746;0.0407031314951731;0.668830467436600;0.000696388699168743];
p=VariableVector(-0.003595083936435,8,{'omega'});
minDist = 0;

%% name states and parameters
xNom=VariableVector(xGue,0,[{'E1r'};{'E1i'};{'n1'};{'E2r'};{'E2i'};{'n2'}]);
alphaNom=VariableVector(param,6,[{'pump1'};{'pump2'}]);

%% define system dynamics
tau=@(x,alpha,p)100/1000;
% lasermod=@(x,xtau,alpha,p)[p(1) * x(2) + 0.5e0 * x(3) * x(1) - 0.20e1 * x(3) * x(2) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(2) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(4) + 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(5);
%                             -p(1) * x(1) + 0.20e1 * x(3) * x(1) + 0.5e0 * x(3) * x(2) - 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(2) - 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(4) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(5);
%                             0.5e-2 * alpha(1) - 0.5e-2 * x(3) - 0.5e-2 * (x(3) + 0.1e1) * (x(1)^2 + x(2)^2);
%                             p(1) * x(5) + 0.5e0 * x(6) * x(4) - 0.20e1 * x(6) * x(5) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(4) + 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(5) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(2);
%                             -p(1) * x(4) + 0.20e1 * x(6) * x(4) + 0.5e0 * x(6) * x(5) - 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(4) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(5) - 0.25e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.25e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(2);
%                             0.5e-2 * alpha(2) - 0.5e-2 * x(6) - 0.5e-2 * (x(6) + 0.1e1) * (x(4)^2 + x(5)^2)  
%                           ];

lasermod=@Laser_DDE_2L;


%% collect them in an object 
aDDE_2L=DDE(@(x,xtau,alpha,p)lasermod(x,xtau,alpha,p),tau,xNom,alphaNom,p);

% clear tau
% clear alphaNom
% clear xGue
% clear param
% clear omega
% clear phi

%% construct object ''aDDENLP''

J=@(x)-abs((x(1)+1i*x(2))+(x(4)+1i*x(5)))^2;

aDDENLP=DDENLP(J,aDDE_2L,xNom,-Inf(6,1),Inf(6,1),[0;0],[0.8;0.5]-sqrt(2)*minDist,-Inf,Inf);
% clear aDDE_1L
% clear J
% clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','active-set','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-14,'TolX',1e-10,'TolCon',1e-10);

%active-set
%% initialize steady state constraints
aDDENLP.initializeStStRot();

%aDDENLP.vars.nominal.x.values

%% initatialize constraints and prepare optimization


aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();



% openfig('Pump_fsolve.fig')
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),aDDENLP.vars.nominal.p.values(1),'ks')
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).p.values(1),'rs')
% 
% figure(1)
% hold on
% plot3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).p.values(1),'rs')
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),aDDENLP.vars.critical(1,1).p.values(1),0.012*aDDENLP.vars.critical(1,1).r.values(1),0.012*aDDENLP.vars.critical(1,1).r.values(2),0)
% 
% figure(1)
% hold on
% quiver3(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),aDDENLP.vars.critical(2,1).p.values(1),0.012*aDDENLP.vars.critical(2,1).r.values(1),0.012*aDDENLP.vars.critical(2,1).r.values(2),0)

openfig('modfold_simple.fig');
xlim([0,1])
ylim([0,1])

figure(1)
hold on
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'ks')

% figure(1)
% hold on
% plot(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),'rs')
% 
% figure(1)
% hold on
% plot(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),'rs')
% 
% figure(1)
% hold on
% quiver(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),0.012*aDDENLP.vars.critical(1,1).r.values(1),0.012*aDDENLP.vars.critical(1,1).r.values(2))
% 
% figure(1)
% hold on
% quiver(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),0.012*aDDENLP.vars.critical(2,1).r.values(1),0.012*aDDENLP.vars.critical(2,1).r.values(2))


%% run optimization
aDDENLP.maxAllowedRealPart=-1;
aDDENLP.minDist=minDist;

aDDENLP.runOptim;
hier = aDDENLP.deconstructOptimum();

figure(1)
hold on
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'kx')

xlabel('A');ylabel('B');
grid on;

%% run simulation

point=aDDENLP.vars.nominal(end);
omegaDist=0*point.p.values(1);

E1=-exp(-0*1i*pi)*(point.x.values(1)+1i*point.x.values(2));
n1=point.x.values(3);
E2=-exp(-0*1i*pi)*(point.x.values(4)+1i*point.x.values(5));
n2=point.x.values(6);


history=@(t)[...
    real(E1*1*exp(+1i*0.5*omegaDist*t));
    imag(E1*1*exp(+1i*0.5*omegaDist*t));
    n1+0*t;
    real(E2*1*exp(-1i*0.5*omegaDist*t));
    imag(E2*1*exp(-1i*0.5*omegaDist*t));
    n2+0*t];



options=ddeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',0.01);
sol=ddesd(aDDENLP,point,history,[0 4],options);

figure(2);clf;
plot(sol.x,sol.y(1:3:end,:))
xlim([0 4]);
 matlab2tikz('2LaserReRotUnStab.tex')

% figure(21);clf;
% 
% tHist=-100:0;
% 
% hist=history(tHist).*exp(1i*omega*tHist);
% Eh1 = (hist(1,:)+1i*hist(2,:)).*exp(-1i*omega*tHist);
% Eh2 = (hist(4,:)+1i*hist(5,:)).*exp(-1i*omega*tHist);
% plot(tHist,real(Eh1),tHist,real(Eh2))
% 


 figure(3);clf;
 A=sol.y(1:3:end,:)+1i*sol.y(2:3:end,:);
 E=A.*repmat(exp(-1i*point.p.values(1)*sol.x*1000),2,1);
 plot(sol.x,real(E))
 xlim([0 4]);
 matlab2tikz('2LaserReFixUnStab.tex')


 
 
 
%% check stability

[~]=aDDENLP.checkStabilityPoint('nominal');
% [~]=aDDENLP.checkStabilityPoint('critical');
% 
% [~]=aDDENLP.checkStabilityAtVertices();
% 
% [~]=aDDENLP.checkStabilityAtRandom();
