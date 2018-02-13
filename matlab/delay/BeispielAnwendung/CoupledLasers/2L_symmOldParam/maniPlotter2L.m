
%% prepare matlab
close all
clear
clc

%% define system to begin with
param=[0.4;0.4];
xGue=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01];

%% name states and parameters
xNom=VariableVector(xGue,0,[{'E1r'};{'E1i'};{'n1'};{'E2r'};{'E2i'};{'n2'}]);
alphaNom=VariableVector(param,6,[{'pump1'};{'pump2'}]);
pNom=VariableVector(-0.02,8,{'omega'});

%% define system dynamics
tau=@(x,alpha,p)100/1000;

lasermod=@Laser_DDE_2L;

%% collect them in an object
aDDE_2L=DDE(@(x,xtau,alpha,p)lasermod(x,xtau,alpha,p),tau,xNom,alphaNom,pNom);

aDDE_2L.modfoldManiHandle = @Mod_Fold_2L;
aDDE_2L.modfoldNVHandle = @NV_Mod_Fold_2L;

% clear tau
% clear alphaNom
% clear xGue
% clear param
% clear omega
% clear phi

%% construct object ''aDDENLP''

minDist=0.01;


J=@(x)-1e2*abs((x(1)+1i*x(2))+(x(4)+1i*x(5))); %Maximale Intensität (mit p1=p2=0.1)

aDDENLP=DDENLP(J,aDDE_2L,xNom,[-Inf;-Inf;-Inf;-Inf;-Inf;-Inf],[Inf;Inf;Inf;Inf;Inf;Inf],[0;0],[1;0.5]-sqrt(2)*minDist,-Inf,Inf);
% clear aDDE_1L
% clear J
% clear xNom

aDDENLP.optionsInitEqCons=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',10000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
aDDENLP.optionsInitOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',200000,'display','off');
aDDENLP.optionsMainOptim=optimoptions('fmincon','Algorithm','interior-point','MaxIter',10000,'MaxFunEvals',200000,'display','iter','TolFun',1e-14,'TolX',1e-10,'TolCon',1e-10);

%active-set
%% initialize steady state constraints
aDDENLP.initializeStSt();

%aDDENLP.vars.nominal.x.values

%% add NV Cons

%Mod-Fold Punkt Branch 1 (100)
xCrit1=VariableVector([0.197742718207656;0.254333015283471;-0.0100548294108757;0.0766499369976581;0.341727571036470;-0.00900159364531876],0,[{'E1rcrit1'};{'E1icrit1'};{'n1crit1'};{'E2rcrit1'};{'E2icrit1'};{'n2crit1'}]);
alphaCrit1=VariableVector([0.106099537864213;0.0870803879148029],0,[{'pump1crit1'};{'pump2crit1'}]);
pCrit1=VariableVector(-0.0193,0,{'omegacrit1'});

aDDENLP.addNVCon('modfold',@Mod_Fold_2L,@NV_Mod_Fold_2L,xCrit1,alphaCrit1,pCrit1);
%aDDENLP.NVCon(end).vars.w1.values=[-2.04600372668366e-06;0.735133515320831;2.69991153805835e-08;0.268291038541200;0.622574199020755;6.34770142102037e-09];
%aus der fsolve-berechnung nur_w_2l
aDDENLP.NVCon(end).vars.w1.values=[-0.202958484289877;0.150220521107475;-0.00334100169657379;-0.505941218404173;0.824774581827560;0.000925615564299419];


% aDDENLP.addNVCon('modfold',@Mod_Fold_2L,@NV_Mod_Fold_2L,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1
clear pCrit1


%Mod-Fold Punkt Branch 1 (100)
xCrit1=VariableVector([0.0766499369976581;0.341727571036470;-0.00900159364531876;0.197742718207656;0.254333015283471;-0.0100548294108757],0,[{'E1rcrit1'};{'E1icrit1'};{'n1crit1'};{'E2rcrit1'};{'E2icrit1'};{'n2crit1'}]);
alphaCrit1=VariableVector([0.0870803879148029;0.106099537864213],0,[{'pump1crit1'};{'pump2crit1'}]);
pCrit1=VariableVector(-0.0193,0,{'omegacrit1'});


aDDENLP.addNVCon('modfold',@Mod_Fold_2L,@NV_Mod_Fold_2L,xCrit1,alphaCrit1,pCrit1);
%aDDENLP.NVCon(end).vars.w1.values=[-2.04600372668366e-06;0.735133515320831;2.69991153805835e-08;0.268291038541200;0.622574199020755;6.34770142102037e-09];
%aus der fsolve-berechnung nur_w_2l
aDDENLP.NVCon(end).vars.w1.values=[-0.505941218404173;0.824774581827560;0.000925615564299419;-0.202958484289877;0.150220521107475;-0.00334100169657379];


% aDDENLP.addNVCon('modfold',@Mod_Fold_2L,@NV_Mod_Fold_2L,xCrit1,alphaCrit1,p);
clear alphaCrit1
clear xCrit1
clear pCrit1



%% initatialize constraints and prepare optimization


aDDENLP.initNVCons();
aDDENLP.concatConstraints();
aDDENLP.concatInitPoints();


figure(1)
hold on
axis equal
plot(aDDENLP.vars.nominal.alpha.values(1),aDDENLP.vars.nominal.alpha.values(2),'ks')

plot(...
    minDist*sqrt(aDDENLP.nAlpha)*sin(0:0.01:2*pi)+aDDENLP.vars.nominal.alpha.values(1),...
    minDist*sqrt(aDDENLP.nAlpha)*cos(0:0.01:2*pi)+aDDENLP.vars.nominal.alpha.values(2),...
    'k')


figure(1)
plot(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),'rs')

figure(1)
plot(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),'rs')
%
figure(1)
quiver(aDDENLP.vars.critical(1,1).alpha.values(1),aDDENLP.vars.critical(1,1).alpha.values(2),0.1*aDDENLP.vars.critical(1,1).r.values(1),0.1*aDDENLP.vars.critical(1,1).r.values(2))

figure(1)
quiver(aDDENLP.vars.critical(2,1).alpha.values(1),aDDENLP.vars.critical(2,1).alpha.values(2),0.1*aDDENLP.vars.critical(2,1).r.values(1),0.1*aDDENLP.vars.critical(2,1).r.values(2))

xlabel('A');ylabel('B');
grid on;



options=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);


for ii = 1:1;%length(aDDENLP.vars.critical)
    x = aDDENLP.vars.critical(ii).x.values;
    alpha = aDDENLP.vars.critical(ii).alpha.values;
    p = aDDENLP.vars.critical(ii).p.values;
    w1 = aDDENLP.vars.critical(ii).w1.values;
    
    currentAlpha = alpha(1);
    
    alphaMani = NaN(2,1);
    alphaMani(:,1) = alpha;
    
    for jj=2:40
        %% prediction
        currentAlpha=currentAlpha - 0.02;
        
        y0 = [x; alpha(2);p;w1];
        
        %% correction
        
        [y,~,exitflag] = fsolve(@(y)aDDENLP.problemDDE.modfoldManiHandle(y(1:6),[currentAlpha,y(7)],y(8),y(9:14)),y0,options);
        
        if exitflag > 0
            x = y(1:6);
            alpha = [currentAlpha;y(7)];
            p = y(8);
            w1 = y(9:14);
        else
            warning('exitflag was %d', exitflag);
            break
        end
        
        alphaMani=[alphaMani,alpha];
    end
    
    x = aDDENLP.vars.critical(ii).x.values;
    alpha = aDDENLP.vars.critical(ii).alpha.values;
    p = aDDENLP.vars.critical(ii).p.values;
    w1 = aDDENLP.vars.critical(ii).w1.values;
    
    currentAlpha = alpha(1);
    
    alphaMani=alphaMani(:,end:-1:1);
    
    
    for jj=1:30
        %% prediction
        currentAlpha=currentAlpha + 0.02;
        
        y0 = [x; alpha(2);p;w1];
        
        %% correction
        
        [y,~,exitflag] = fsolve(@(y)aDDENLP.problemDDE.modfoldManiHandle(y(1:6),[currentAlpha,y(7)],y(8),y(9:14)),y0,options);
        
        if exitflag > 0
            x = y(1:6);
            alpha = [currentAlpha;y(7)];
            p = y(8);
            w1 = y(9:14);
        else
            warning('exitflag was %d', exitflag);
            alphaMani=[alphaMani,alpha];
            break
        end
        alphaMani=[alphaMani,alpha];
    end
    
    
end

plot([alphaMani(1,1),alphaMani(1,end)],[alphaMani(2,1),alphaMani(2,end)],'m');
plot([alphaMani(2,1),alphaMani(2,end)],[alphaMani(1,1),alphaMani(1,end)],'m');




