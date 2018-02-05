% % Berechung des Verlaufes der Kostenfunktion im Beispiel drei asymmetrisch
% % gekoppelter Laser --> Abbildung in 3D
% 
% clc;
% clear;
% close all;
% Sol1=[];
% Sol2=[];
% Sol3=[];
% 
% var0=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01];
% omega=-0.02;
% 
% p1=40;
% for p2=2:2:40
%     for p3=2:2:40
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_Setup1_fkt(var,pump,omega);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)))^2; 
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol1=[Sol1;Erg];
%     end
% end
% 
% 
% for p1=2:2:40
%     p2=40;
%     for p3=2:2:40
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_Setup1_fkt(var,pump,omega);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)))^2;
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol2=[Sol2;Erg];
%     end
% end
% 
% 
% for p1=2:2:40
%     for p2=2:2:40
%         p3=40;
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_Setup1_fkt(var,pump,omega);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)))^2; 
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol3=[Sol3;Erg];
%         
%     end
% end
% 
% 
% figure(1)%;clf;
% hold on
% scatter3(Sol1(:,1),Sol1(:,2),Sol1(:,3),150,Sol1(:,4),'filled')
% scatter3(Sol2(:,1),Sol2(:,2),Sol2(:,3),150,Sol2(:,4),'filled')
% scatter3(Sol3(:,1),Sol3(:,2),Sol3(:,3),150,Sol3(:,4),'filled')
% colorbar
% xlabel('A')
% ylabel('B')
% zlabel('C')
% axis equal
% grid on
% xlim([0 0.5])
% ylim([0 0.5])
% zlim([0 0.5])
% 


% Berechnung des Verlaufes der Kostenfunktion im System von drei
% symmetrisch gekoppelten Lasern

clc;
clear;
close all;

var0=[0.1910;0.2732;-0.01;0.1910;0.2732;-0.01;0.1910;0.2732;-0.01];
omega=-0.02;

%% in p1-p2-plane
steplength=0.005;

p1=0:steplength:1;
p2=0:steplength:1;
p3 = 0.3827;

cost = NaN(length(p1),length(p2));


for ii=1:length(p1)
    for jj=1:length(p2)
        pump = [p1(ii) p2(jj) p3];
        f =@(var) DDE_3L_SETUP_1(var(1:9),var(1:9),pump,var(end));
        opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
        [out]=fsolve(f,[var0;omega],opts);
        Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)));
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol1=[Sol1;Erg];
        cost(ii,jj)=Guete;
    end
end

% for p1=1:2:40
%     p2=40;
%     for p3=1:2:40
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_fkt(var,pump);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)));
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol2=[Sol2;Erg];
%     end
% end
% 
% for p1=1:2:40
%     for p2=1:2:40
%         p3=40;
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_fkt(var,pump);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)));
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol3=[Sol3;Erg];
%     end
% end


figure(2);clf;
hold on
axis equal
box on
contour(p1,p2,cost',0:0.1:3)
xlim([0.2 0.9])
ylim([0.2 0.9])
grid on
xlabel('p_1');
ylabel('p_2');

matlab2tikz('contourp1p2asym3L');

return

%% in p2-p3-plane
steplength=0.01;

p1 = 0.5177;
p2 = 0:steplength:1;
p3 = 0:steplength:1;

cost = NaN(length(p2),length(p3));


for ii=1:length(p2)
    for jj=1:length(p3)
        pump = [p1 p2(ii) p3(jj)];
        f =@(var) DDE_3L_SETUP_1(var(1:9),var(1:9),pump,var(end));
        opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
        [out]=fsolve(f,[var0;omega],opts);
        Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)));
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol1=[Sol1;Erg];
        cost(ii,jj)=Guete;
    end
end

% for p1=1:2:40
%     p2=40;
%     for p3=1:2:40
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_fkt(var,pump);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)));
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol2=[Sol2;Erg];
%     end
% end
% 
% for p1=1:2:40
%     for p2=1:2:40
%         p3=40;
%         pump = [p1*0.01 p2*0.01 p3*0.01];
%         f =@(var) Syncmanifold_3L_fkt(var,pump);
%         opts=optimoptions('fsolve','Algorithm','levenberg-marquardt');
%         [out]=fsolve(f,var0,opts);
%         Guete=abs((out(1)+1i*out(2))+(out(4)+1i*out(5))+(out(7)+1i*out(8)));
%         Erg=[pump(1),pump(2),pump(3),Guete];
%         Sol3=[Sol3;Erg];
%     end
% end


figure(2);clf;
hold on
axis equal
box on
contour(p2,p3,cost',0:0.05:3)
xlim([0 0.8])
ylim([0 0.8])
grid on
xlabel('p_2');
ylabel('p_3');

matlab2tikz('contourp2p3asym3L');


