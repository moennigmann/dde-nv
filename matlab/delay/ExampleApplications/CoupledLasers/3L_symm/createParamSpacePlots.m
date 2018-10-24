close all
clear
clc

load('opt3symLaser.mat')
alphaNom=[0.517713706146493;0.482679491924311;0.482679491924311];

%% plot p1-p2 plane
figure(1);clf;
hold on
box on
grid on
axis equal

plot([0 1],[0.5,0.5],'k')
plot([0.8 0.8],[0,1],'k')
plot(alphaNom(1),alphaNom(2),'x')

xlim([0 1]);
ylim([0 1]);

for ii=1:3
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,~]=plane_intersect(r,alphaCrit,[0 0 1]',[0 0 alphaNom(3)]);
    n=n*sign(n(2));
    plot([p(1),p(1)+10*n(1)],[p(2),p(2)+10*n(2)],'k')
end

% plot([ 0.5277 0.5277 0.5077 0.5077 0.5277 ],[0.4727 0.4927 0.4927 0.4727 0.4727],'b')
% plot(sqrt(3)*0.01*cos(2*pi*(0:0.01:1))+0.5177,sqrt(3)*0.01*sin(2*pi*(0:0.01:1))+0.4827,'r')

xlabel('p_1')
ylabel('p_2')


% matlab2tikz('3LsymParamP1P2.tex')

%% plot p2-p3 plane
figure(2);clf;
hold on
box on
grid on
axis equal

plot([0 1],[0.5,0.5],'k')
plot([0.5,0.5],[0 1],'k')
plot(alphaNom(2),alphaNom(3),'x')

xlim([0 1]);
ylim([0 1]);

for ii=1:1
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,check]=plane_intersect(r,alphaCrit,[1 0 0]',[alphaNom(1) 0 0]);
    n=n*sign(n(3));
    
    plot([p(2),p(2)+10*n(2)],[p(3),p(3)+10*n(3)],'k')
end

xlabel('p_2')
ylabel('p_3')

matlab2tikz('3LsymParamP2P3.tex')