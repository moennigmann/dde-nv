% Plotten der kritischen Mannigfaltigkeiten (Berechnung nicht enthalten)
% 3L-asymmetrisch
%
% Jens Müller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

%Hier kann entschieden werden, ob die Lösungszweige mitgeplottet werden
%oder nicht

%openfig('Stab_3L_fixed_K.fig')

figure(1)
hold on
grid on

%% Berechnung der Schnittvektoren

%Normalenvektoren, bekannt aus der Optimierung
NV1=[0.727532830026351,3.65252886824957e-21,-0.686072868749267];
%Stützvektoren aus der numerischen Berchnung (wähle Punkt 10)
P1=[0.0552028416137532;0.0470075845635549;0.0596916770582316];

NV2=[0.727532830026354,-0.686072868749264,-1.09521933437694e-21];
P2=[0.0510141762018182;0.0552509135326719;0.0550960738310083];

NV3=[-0.0275066240363199;0.0421637255044102;5.88293483917290e-10];
P3=[0.584998675541705;0.357453443197026;0.506111304487198];

NV4=[-0.0124117167910764,4.62019480599190e-33,0.0217073850691745];
P4=[0.584998675537902;0.506111304496207;0.357453443194848];

%Testwerte
%aNV3=[-0.0124117167910057;0.0217073850690769;-5.23617105121988e-27];
%P3=[0.0675407567864075;0.0315899466133940;0.0715490447761822];
%NV4=[-0.0133146582293309,-3.90167472311009e-26,0.0229520992832143];
%P4=[0.0670752824844416;0.0479271839661028;0.0312329568134677];

%Plotten der Ebenen zur Überprüfung

%      w = null(NV1); % Find two orthonormal vectors which are orthogonal to v
%    [P,Q] = meshgrid(0:5); % Provide a gridwork (you choose the size)
%    X = P1(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
%    Y = P1(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
%    Z = P1(3)+w(3,1)*P+w(3,2)*Q;
%    surf(X,Y,Z)
%       xlim([0 5])
%     ylim([0 5])
%      zlim([0 5])
%   hold on
%
%        w = null(NV2); % Find two orthonormal vectors which are orthogonal to v
%    [P,Q] = meshgrid(0:5); % Provide a gridwork (you choose the size)
%    X = P2(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
%    Y = P2(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
%    Z = P2(3)+w(3,1)*P+w(3,2)*Q;
%    surf(X,Y,Z)
%       xlim([0 5])
%     ylim([0 5])
%      zlim([0 5])
%
%           w = null(NV3); % Find two orthonormal vectors which are orthogonal to v
%    [P,Q] = meshgrid(0:5); % Provide a gridwork (you choose the size)
%    X = P3(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
%    Y = P3(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
%    Z = P3(3)+w(3,1)*P+w(3,2)*Q;
%    surf(X,Y,Z)
%       xlim([0 5])
%     ylim([0 5])
%      zlim([0 5])
%
%           w = null(NV4); % Find two orthonormal vectors which are orthogonal to v
%    [P,Q] = meshgrid(0:5); % Provide a gridwork (you choose the size)
%    X = P4(1)+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
%    Y = P4(2)+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
%    Z = P4(3)+w(3,1)*P+w(3,2)*Q;
%    surf(X,Y,Z)
%
%    xlim([0 5])
%     ylim([0 5])
%      zlim([0 5])


%Schnittgeraden berechnetn mit Plane-Intersect
[GP1,GN1,check1]=plane_intersect(NV1,P1,NV2,P2);

[GP2,GN2,check2]=plane_intersect(NV1,P1,NV3,P3);

[GP3,GN3,check3]=plane_intersect(NV3,P3,NV4,P4);

[GP4,GN4,check4]=plane_intersect(NV4,P4,NV2,P2);


%Normieren der Richtungsvektoren

nGN1=(GN1/norm(GN1))*5;
nGN2=(GN2/norm(GN2))*5;
nGN3=(GN3/norm(GN3))*5;
nGN4=(GN4/norm(GN4))*5;


% S2n=(S2/norm(S2))/1.6;
% S3n=(S3/norm(S3))/1.6;
% S4n=(S4/norm(S4))/1.6;

% quiver3(0,0,0,-S1n(1),-S1n(2),-S1n(3));
% quiver3(0,0,0,S2n(1),S2n(2),S2n(3));
% quiver3(0,0,0,S3n(1),S3n(2),S3n(3));
% quiver3(0,0,0,S4n(1),S4n(2),S4n(3));

%Plotten der Schnittgeraden
%
% plot3([GP1(1);-GP1(1)-nGN1(1)],[GP1(2);-GP1(2)-nGN1(2)],[GP1(3);-GP1(3)-nGN1(3)],'k','LineWidth',1);
% plot3([GP2(1);+GP2(1)+nGN2(1)],[GP2(2);+GP2(2)+nGN2(2)],[GP2(3);+GP2(3)+nGN2(3)],'k','LineWidth',1);
% plot3([GP3(1);+GP3(1)+nGN3(1)],[GP3(2);+GP3(2)+nGN3(2)],[GP3(3);+GP3(3)+nGN3(3)],'k','LineWidth',1);
% plot3([GP4(1);+GP4(1)+nGN4(1)],[GP4(2);+GP4(2)+nGN4(2)],[GP4(3);+GP4(3)+nGN4(3)],'k','LineWidth',1);

%Plotten der Schnittgeraden getrickst

plot3([0;-GP1(1)-nGN1(1)],[0;-GP1(2)-nGN1(2)],[0;-GP1(3)-nGN1(3)],'k','LineWidth',1);
plot3([0;+GP2(1)+nGN2(1)],[0;+GP2(2)+nGN2(2)],[0;+GP2(3)+nGN2(3)],'k','LineWidth',1);
plot3([0;+GP3(1)+nGN3(1)],[0;+GP3(2)+nGN3(2)],[0;+GP3(3)+nGN3(3)],'k','LineWidth',1);
plot3([0;+GP4(1)+nGN4(1)],[0;+GP4(2)+nGN4(2)],[0;+GP4(3)+nGN4(3)],'k','LineWidth',1);


% %Plotten der transparenten Mannigfaltigkeiten
% 
% x1=[GP1(1)  -GP1(1)-nGN1(1) GP2(1)+nGN2(1) GP2(1)]
% y1=[GP1(2)  -GP1(2)-nGN1(2) GP2(2)+nGN2(2) GP2(2)]
% z1=[GP1(3)  -GP1(3)-nGN1(3) GP2(3)+nGN2(3) GP2(3)]
% patch(x1,y1,z1,'r','LineStyle','none','FaceAlpha',0.5);
% 
% x2=[GP1(1)  -GP1(1)-nGN1(1) GP4(1)+nGN4(1) GP4(1)]
% y2=[GP1(2)  -GP1(2)-nGN1(2) GP4(2)+nGN4(2) GP4(2)]
% z2=[GP1(3)  -GP1(3)-nGN1(3) GP4(3)+nGN4(3) GP4(3)]
% patch(x2,y2,z2,'r','LineStyle','none','FaceAlpha',0.5);
% 
% x3=[GP2(1)  GP2(1)+nGN2(1) GP3(1)+nGN3(1) GP3(1)]
% y3=[GP2(2)  GP2(2)+nGN2(2) GP3(2)+nGN3(2) GP3(2)]
% z3=[GP2(3)  GP2(3)+nGN2(3) GP3(3)+nGN3(3) GP3(3)]
% patch(x3,y3,z3,'b','LineStyle','none','FaceAlpha',0.5);
% 
% x4=[GP4(1)  GP4(1)+nGN4(1) GP3(1)+nGN3(1) GP3(1)]
% y4=[GP4(2)  GP4(2)+nGN4(2) GP3(2)+nGN3(2) GP3(2)]
% z4=[GP4(3)  GP4(3)+nGN4(3) GP3(3)+nGN3(3) GP3(3)]
% patch(x4,y4,z4,'b','LineStyle','none','FaceAlpha',0.5);

%Plotten der transparenten Mannigfaltigkeiten
% getrickst
x1=[0 -GP1(1)-nGN1(1) GP2(1)+nGN2(1)];
y1=[0 -GP1(2)-nGN1(2) GP2(2)+nGN2(2)];
z1=[0  -GP1(3)-nGN1(3) GP2(3)+nGN2(3)];
patch(x1,y1,z1,'r','LineStyle','none','FaceAlpha',0.5);

x2=[0 -GP1(1)-nGN1(1) GP4(1)+nGN4(1)];
y2=[0 -GP1(2)-nGN1(2) GP4(2)+nGN4(2)];
z2=[0 -GP1(3)-nGN1(3) GP4(3)+nGN4(3)];
patch(x2,y2,z2,'r','LineStyle','none','FaceAlpha',0.5);

x3=[0 GP2(1)+nGN2(1) GP3(1)+nGN3(1)];
y3=[0 GP2(2)+nGN2(2) GP3(2)+nGN3(2)];
z3=[0 GP2(3)+nGN2(3) GP3(3)+nGN3(3)];
patch(x3,y3,z3,'b','LineStyle','none','FaceAlpha',0.5);

x4=[0 GP4(1)+nGN4(1) GP3(1)+nGN3(1)];
y4=[0 GP4(2)+nGN4(2) GP3(2)+nGN3(2)];
z4=[0 GP4(3)+nGN4(3) GP3(3)+nGN3(3)];
patch(x4,y4,z4,'b','LineStyle','none','FaceAlpha',0.5);

axis equal
box on
xlim([0 0.9])
ylim([0 0.9])
zlim([0 0.9])
xlabel('p1')
ylabel('p2')
zlabel('p3')
view([-30 30])
% matlab2tikz('3LasymManis.tex')



%% load data for cuts
load('opt3asymLaser.mat')

alphaNom=aDDENLP.vars.nominal.alpha.values;
plot3(alphaNom(1), alphaNom(2), alphaNom(3), 'kx')

load('init3asymLaser.mat','aDDENLP')

%% plot p1-p2 plane cuts


for ii=1:4
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,~]=plane_intersect(r,alphaCrit,[0 0 1]',[0 0 alphaNom(3)]);
    
    hp3(ii)=plot3([p(1)-10*n(1),p(1)+10*n(1)],[p(2)-10*n(2),p(2)+10*n(2)],[alphaNom(3),alphaNom(3)])
end

% matlab2tikz('3LasymParamP1P2Manis.tex')

%% plot p2-p3 plane cuts

for ii=1:4
    alphaCrit=aDDENLP.vars.critical(ii).alpha.values;
    r=aDDENLP.vars.critical(ii).r.values;
    r=r/norm(r,2);
    
    [p,n,~]=plane_intersect(r,alphaCrit,[1 0 0]',[alphaNom(1) 0 0 ]);
    
    hp1(ii)=plot3([alphaNom(1),alphaNom(1)],[p(2)-10*n(2),p(2)+10*n(2)],[p(3)-10*n(3),p(3)+10*n(3)])
end

% matlab2tikz('3LasymParamP2P3Manis.tex')

