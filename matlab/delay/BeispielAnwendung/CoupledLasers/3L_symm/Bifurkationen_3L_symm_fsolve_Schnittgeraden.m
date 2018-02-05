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
  NV1=[-0.785483832000789;0.437615784487579;0.437615784487578];
%Stützvektoren aus kritischen Punkten
  P1=[0.439197280946214;0.392420272270734;0.392420272270734];
  
  NV2=[0.445558320355562;-0.785391553005577;0.429695114737671];
  P2=[0.403818644331403;0.438380714586315;0.379001549615898];
    
  NV3=[0.445593631346356;0.429657770418340;-0.785391950571629];
  P3=[0.403816937463878;0.379003392617136;0.438380701021090];


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

[GP2,GN2,check2]=plane_intersect(NV2,P2,NV3,P3);

[GP3,GN3,check3]=plane_intersect(NV3,P3,NV1,P1);



%Normieren der Richtungsvektoren

nGN1=(GN1/norm(GN1))*2;
nGN2=(GN2/norm(GN2))*2;
nGN3=(GN3/norm(GN3))*2;



% S2n=(S2/norm(S2))/1.6;
% S3n=(S3/norm(S3))/1.6;
% S4n=(S4/norm(S4))/1.6;

% quiver3(0,0,0,-S1n(1),-S1n(2),-S1n(3));
% quiver3(0,0,0,S2n(1),S2n(2),S2n(3));
% quiver3(0,0,0,S3n(1),S3n(2),S3n(3));
% quiver3(0,0,0,S4n(1),S4n(2),S4n(3));

%Plotten der Schnittgeraden
% 
plot3([GP1(1);+GP1(1)+nGN1(1)],[GP1(2);+GP1(2)+nGN1(2)],[GP1(3);+GP1(3)+nGN1(3)],'k','LineWidth',1);
plot3([GP2(1);+GP2(1)+nGN2(1)],[GP2(2);+GP2(2)+nGN2(2)],[GP2(3);+GP2(3)+nGN2(3)],'k','LineWidth',1);
plot3([GP3(1);+GP3(1)+nGN3(1)],[GP3(2);+GP3(2)+nGN3(2)],[GP3(3);+GP3(3)+nGN3(3)],'k','LineWidth',1);

%Plotten der Schnittgeraden getrickst

% plot3([0;-GP1(1)-nGN1(1)],[0;-GP1(2)-nGN1(2)],[0;-GP1(3)-nGN1(3)],'k','LineWidth',1);
% plot3([0;+GP2(1)+nGN2(1)],[0;+GP2(2)+nGN2(2)],[0;+GP2(3)+nGN2(3)],'k','LineWidth',1);
% plot3([0;+GP3(1)+nGN3(1)],[0;+GP3(2)+nGN3(2)],[0;+GP3(3)+nGN3(3)],'k','LineWidth',1);


%Plotten der transparenten Mannigfaltigkeiten

x1=[GP1(1)  GP1(1)+nGN1(1) GP3(1)+nGN3(1) GP3(1)];
y1=[GP1(2)  GP1(2)+nGN1(2) GP3(2)+nGN3(2) GP3(2)];
z1=[GP1(3)  GP1(3)+nGN1(3) GP3(3)+nGN3(3) GP3(3)];
patch(x1,y1,z1,'r','LineStyle','none','FaceAlpha',0.5);

x2=[GP1(1)  GP1(1)+nGN1(1) GP2(1)+nGN2(1) GP2(1)];
y2=[GP1(2)  GP1(2)+nGN1(2) GP2(2)+nGN2(2) GP2(2)];
z2=[GP1(3)  GP1(3)+nGN1(3) GP2(3)+nGN2(3) GP2(3)];
patch(x2,y2,z2,'r','LineStyle','none','FaceAlpha',0.5);

x3=[GP2(1)  GP2(1)+nGN2(1) GP3(1)+nGN3(1) GP3(1)];
y3=[GP2(2)  GP2(2)+nGN2(2) GP3(2)+nGN3(2) GP3(2)];
z3=[GP2(3)  GP2(3)+nGN2(3) GP3(3)+nGN3(3) GP3(3)];
patch(x3,y3,z3,'r','LineStyle','none','FaceAlpha',0.5);


% %Plotten der transparenten Mannigfaltigkeiten
% % getrickst
% x1=[0 -GP1(1)-nGN1(1) GP2(1)+nGN2(1)]
% y1=[0 -GP1(2)-nGN1(2) GP2(2)+nGN2(2)]
% z1=[0  -GP1(3)-nGN1(3) GP2(3)+nGN2(3)]
% patch(x1,y1,z1,'r','LineStyle','none','FaceAlpha',0.5);
% 
% x2=[0 -GP1(1)-nGN1(1) GP4(1)+nGN4(1)]
% y2=[0 -GP1(2)-nGN1(2) GP4(2)+nGN4(2)]
% z2=[0 -GP1(3)-nGN1(3) GP4(3)+nGN4(3)]
% patch(x2,y2,z2,'r','LineStyle','none','FaceAlpha',0.5);
% 
% x3=[0 GP2(1)+nGN2(1) GP3(1)+nGN3(1)]
% y3=[0 GP2(2)+nGN2(2) GP3(2)+nGN3(2)]
% z3=[0 GP2(3)+nGN2(3) GP3(3)+nGN3(3)]
% patch(x3,y3,z3,'b','LineStyle','none','FaceAlpha',0.5);
% 
% x4=[0 GP4(1)+nGN4(1) GP3(1)+nGN3(1)]
% y4=[0 GP4(2)+nGN4(2) GP3(2)+nGN3(2)]
% z4=[0 GP4(3)+nGN4(3) GP3(3)+nGN3(3)]
% patch(x4,y4,z4,'b','LineStyle','none','FaceAlpha',0.5);

axis equal
box on
xlim([0 0.5178])
ylim([0 0.9])
zlim([0 0.9])
xlabel('p1')
ylabel('p2')
zlabel('p3')

view([-30 30])

plot3(0.5177,0.4827,0.4827,'kx');

patch([0 0.9 0.9 0],[0 0 0.9 0.9],[0.5 0.5 0.5 0.5],'k','FaceAlpha',0.1)
patch([0 0.9 0.9 0],[0.5 0.5 0.5 0.5],[0 0 0.9 0.9],'k','FaceAlpha',0.1)
patch([0.8 0.8 0.8 0.8],[0 0.9 0.9 0],[0 0 0.9 0.9],'k','FaceAlpha',0.1)

patch([0.5177 0.5177 0.5177 0.5177],[0 0.9 0.9 0],[0 0 0.9 0.9],'m','FaceAlpha',0.1)
patch([0 0.9 0.9 0],[0 0 0.9 0.9],[0.4827 0.4827 0.4827 0.4827],'m','FaceAlpha',0.1)




% matlab2tikz('3LsymManis.tex')
