%% This script combines all CSTR Optimizations to compare them afterwards

%% clean up
close all
clear all
clc

runOptimFlag=1;

%% run the scripts
if runOptimFlag

disp('starting with first system: one CSTR')
tic
oneCSTROptim
optimTime1=toc;

save('allOptimResultsNCSTR.mat','oneCSTRDDENLP')
save('allOptimResultsNCSTR.mat','optimTime1','-append')

disp('starting with second system: two CSTRs')
tic
twoCSTROptim
optimTime2=toc;

save('allOptimResultsNCSTR.mat','twoCSTRDDENLP','-append')
save('allOptimResultsNCSTR.mat','optimTime2','-append')

disp('starting with third system: three CSTRs')
tic
threeCSTROptim
optimTime3=toc;

save('allOptimResultsNCSTR.mat','threeCSTRDDENLP','-append')
save('allOptimResultsNCSTR.mat','optimTime3','-append')

close all
clear all
end


%% load results for further treatment

load('allOptimResultsNCSTR.mat')


%% run simulations

% 1 CSTR
figure(1);clf;
point=oneCSTRDDENLP.vars.nominal(1);

history=abs(point.x.values).*[0.9;1.1;1;1;1;1];
options=ddeset();
sol=ddesd(oneCSTRDDENLP,point,history,[0 100],options);

hold on
plot(sol.x,sol.y(1:3:end,:))
plot(sol.x,sol.y(2:3:end,:),'-.')
plot([sol.x(1,1),sol.x(1,end)],[point.x.values(1:3:end),point.x.values(1:3:end)]','--')

% 2 CSTR
figure(2);clf;
point=twoCSTRDDENLP.vars.nominal(1);

history=abs(point.x.values).*[0.9;1.1;1;1;1;1;1;1;1];
options=ddeset();
sol=ddesd(twoCSTRDDENLP,point,history,[0 100],options);

hold on
plot(sol.x,sol.y(1:3:end,:))
plot(sol.x,sol.y(2:3:end,:),'-.')
plot([sol.x(1,1),sol.x(1,end)],[point.x.values(1:3:end),point.x.values(1:3:end)]','--')


% 3 CSTR
figure(3);clf;
point=threeCSTRDDENLP.vars.nominal(1);

history=abs(point.x.values).*[0.9;1.1;1;1;1;1;1;1;1;1;1;1];
options=ddeset();
sol=ddesd(threeCSTRDDENLP,point,history,[0 100],options);

hold on
plot(sol.x,sol.y(1:3:end,:))
plot(sol.x,sol.y(2:3:end,:),'-.')
plot([sol.x(1,1),sol.x(1,end)],[point.x.values(1:3:end),point.x.values(1:3:end)]','--')



%% comparison of results
Jopt(1) = oneCSTRDDENLP.optJ;
Jopt(2) = twoCSTRDDENLP.optJ;
Jopt(3) = threeCSTRDDENLP.optJ;

figure(123);clf;
bar([1,2,3],-1/100*[Jopt(1),Jopt(2),Jopt(3)])
hold on
plot([1 3],-1/100*[Jopt(1),Jopt(3)],'ro-')
plot([2],-1/100*[Jopt(2)],'ro-')