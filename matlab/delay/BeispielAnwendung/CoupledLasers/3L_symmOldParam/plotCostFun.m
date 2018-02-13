%% creates a contour plot of the cost function


%% clean up
close all
clear
clc

recalc=1;

if recalc==1 || ~exist('plotCostFun3Lsym.mat','file')
    
    %% define Cost function
    
    J = @(x)-abs((x(1)+1i*x(2))+(x(4)+1i*x(5)+(x(7)+1i*x(8))))^2;
    
    
    %% define range
    step = 0.1;
    
    p1 = 0:step:0.9;
    p2 = 0:step:0.6;
    p3 = 0.482679491924311;
    
    
    %% eval cost function
    
    
    cost = NaN(length(p1),length(p2));
    omega = NaN(size(cost));
    xCol = NaN(length(p1),length(p2),10);
    
    
    x01=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;-0.0];
    x02=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;-0.02];
    x01=reshape(x01,1,1,10);
    x02=reshape(x02,1,1,10);
    
    
    options=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',200000,'display','off','TolFun',1e-7,'TolX',1e-7);
    
    
    for ii=length(p1):-1:1
        for jj = length(p2):-1:1
            
            rhs = @(x)DDE_3L_symm(x(1:9),x(1:9),[p1(ii);p2(jj);p3],x(7));
            [x1,~,exitflag1] = fsolve(rhs,x01,options);
            [x2,~,exitflag2] = fsolve(rhs,x02,options);
            
            if J(x1)<J(x2)
                x=x1;
            else
                x=x2;
            end
            
            
            if exitflag1 > 0 || exitflag2 > 0
                cost(ii,jj) = J(x);
                omega(ii,jj) = x(1,1,7);
                xCol(ii,jj,:)=x;
            else
                warning('not successful,exitflags were %d and %d',exitflag1,exitflag2)
            end
        end
    end
    
    save('plotCostFun3Lsym.mat')
    
    
else
    
    load('plotCostFun3Lsym.mat')
end

%% create plot
figure(1);
clf;
hold on
% axis equal
grid on
box on

xlim([0 p1(end)])
ylim([0 p2(end)])


% countour(p1,p2,-cost',0:0.1:3);
surf(p1,p2,-cost');
% axis equal

% matlab2tikz('cost2L.tex')

% figure(2);
% clf;
% hold on
% % axis equal
% grid on
% box on
% 
% xlim([0 p1(end)])
% ylim([0 p2(end)])
% 
% surf(p1,p2,-cost');
% view(3)



