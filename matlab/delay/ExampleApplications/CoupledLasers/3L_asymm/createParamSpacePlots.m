close all
clear
clc

load('opt3asymLaser.mat')

alphaNom = [0.594273236532978;0.606081751647561;0.382679491924311];

%% plot p1-p2 plane
figure(1);clf;
hold on
box on
grid on
axis equal

plot([0 1],[0.7,0.7],'k')
plot([0.8 0.8],[0,1],'k')
plot(alphaNom(1),alphaNom(2),'x')

xlim([0 1.5]);
ylim([0 1]);

for ii=1:2
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,~]=plane_intersect(r,alphaCrit,[0 0 1]',[0 0 alphaNom(3)]);
    
    plot([p(1)-10*n(1),p(1)+10*n(1)],[p(2)-10*n(2),p(2)+10*n(2)])
end

for ii=3:4
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,~]=plane_intersect(r,alphaCrit,[0 0 1]',[0 0 alphaNom(3)]);
    
    plot([p(1)-10*n(1),p(1)+10*n(1)],[p(2)-10*n(2),p(2)+10*n(2)],'--')
end


xlabel('p_1')
ylabel('p_2')

% matlab2tikz('3LasymParamP1P2.tex')


%% plot p2-p3 plane
figure(2);clf;
hold on
box on
grid on
axis equal

plot([0 1],[0.4,0.4],'k')
plot([0.7,0.7],[0 1],'k')
plot(alphaNom(2),alphaNom(3),'x')

xlim([0 1]);
ylim([0 1]);

for ii=1:2
    
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,check]=plane_intersect(r,alphaCrit,[1 0 0]',[alphaNom(1) 0 0]);
    
    plot([p(2)-10*n(2),p(2)+10*n(2)],[p(3)-10*n(3),p(3)+10*n(3)])
    
end
for ii=3:4
    
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,check]=plane_intersect(r,alphaCrit,[1 0 0]',[alphaNom(1) 0 0]);
    
    plot([p(2)-10*n(2),p(2)+10*n(2)],[p(3)-10*n(3),p(3)+10*n(3)],'--')
    
end

xlabel('p_2')
ylabel('p_3')

% matlab2tikz('3LasymParamP2P3.tex')