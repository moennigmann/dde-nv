% In diesem Skript werden die Bifurkationen für 3 asymmetrisch gekoppelte
% Laser gezeichnet. Dazu werden numerische Ergebnisse auch vorherigen
% Berechnungen verwendet.
%
% 22.02.2018 Jens Müller

clc;
clear;
close all;

% Diese Matrix enhält die Berechneten Werte der Krümmung der mod. Hopf
A=[0.124214322239783 - 1.74636492082703e-25i;0.154092433509132 - 4.37007765443663e-28i;0.183005236384061 - 2.44173160010315e-28i;0.211019943660515 - 1.18509958771408e-28i;0.238195334905376 - 4.89849806933401e-29i;0.264584095647409 - 3.71078543781371e-28i;0.290233980593492 - 9.05688300143055e-26i;0.315188538292515 - 1.97112938868311e-24i;0.339487635660642 - 1.03348241397579e-24i;0.363167873263014 - 7.29435973558964e-24i;0.386262931185838 + 8.51673574560843e-26i;0.408803865375825 + 2.87376315240620e-25i;0.430819365687786 + 1.70316243132177e-24i;0.452335982774569 - 2.01924998899196e-24i;0.473378328813576 - 6.27232148066104e-24i;0.493969255836324 + 3.20100918843940e-23i;0.514130014655693 - 1.98395319929747e-24i];

% Dieser Normalenvektor beschreibt die Schnittgerade zwischen den beiden
% mod. Fold
S1n=[-0.554783237083809;-0.588309255345697;-0.588309255345692]
% und seine Steigung
m=S1n(2)/S1n(1)

% Ablegen der nötigen Daten in Variablen
q=[0.2:0.05:1]';
w=q*m;
p=A(:,1);

figure(1)

hold on

plot3(q,p,p,'k','LineWidth',2)

plot3(q,w,p,'k','LineWidth',2)

plot3(q,w,w,'k','LineWidth',2)

plot3(q,p,w,'k','LineWidth',2)

%% Get patches
% first fold patches
for i=1:1:16
    patch([q(i), q(i+1), q(i+1), q(i)],[p(i), p(i+1), w(i+1), w(i)],[w(i), w(i+1), w(i+1), w(i)],'r','LineStyle','none','FaceAlpha',0.5);
end
% second fold patches
for i=1:1:16
    patch([q(i), q(i+1), q(i+1), q(i)],[w(i), w(i+1), w(i+1), w(i)],[w(i), w(i+1), p(i+1), p(i)],'r','LineStyle','none','FaceAlpha',0.5);
end
% first hopf patches
for i=1:1:16
    patch([q(i), q(i+1), q(i+1), q(i)],[w(i), w(i+1), p(i+1), p(i)],[p(i), p(i+1), p(i+1), p(i)],'b','LineStyle','none','FaceAlpha',0.5);
end
% second hopf patches
for i=1:1:16
    patch([q(i), q(i+1), q(i+1), q(i)],[p(i), p(i+1), p(i+1), p(i)],[w(i), w(i+1), p(i+1), p(i)],'b','LineStyle','none','FaceAlpha',0.5);
end



