function [ rhs ] = supplyChainModel(x,xtau,alpha,~)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

xx=[x,xtau];

%% parameters
% getSupplyChainParameterVector;

T1    = 1;
aWIP1 = 0.01;
ai1   = 0.1;
hD1   = alpha(1);
hP1   = alpha(2);
hA1   = 11;

T2    = 1;
aWIP2 = 0.01;
ai2   = 0.1;
hD2   = alpha(3);
hP2   = alpha(4);
hA2   = 12;

T3    = 1;
aWIP3 = 0.01;
ai3   = 0.1;
hD3   = alpha(5);
hP3   = alpha(6);
hA3   = 13;


uu=ones(1,13)*alpha(7);



%% Matrices
A0 = zeros(6);
A0(1,2)=1;
A0(2,2)=-1/T1;
A0(3,4)=1;
A0(4,4)=-1/T2;
A0(5,6)=1;
A0(6,6)=-1/T3;


AD1 = zeros(6);
AD1(2,1) = -aWIP1/T1;


AP1 = zeros(6);
AP1(2,1) = aWIP1/T1;


AT1 = zeros(6);
AT1(2,1) = -ai1/T1;

AA1=zeros(6);


AD2 = zeros(6);
AD2(4,1) = (ai2*hA2 - (1+aWIP2*hP2))/(T2*hA2); 
AD2(4,3) = -aWIP2/T2;


AP2 = zeros(6);
AP2(4,3) = aWIP2/T2;


AT2 = zeros(6);
AT2(4,3) = -ai2/T2;


AA2 = zeros(6);
AA2(4,1) = (1+aWIP2*hP2)/(T2*hA2);


AD3 = zeros(6);
AD3(6,3) = (ai3*hA3 - (1+aWIP3*hP3))/(T3*hA3);
AD3(6,5) = -aWIP3/T3;


AP3 = zeros(6);
AP3(6,5) = aWIP3/T3;


AT3 = zeros(6);
AT3(6,5) = -ai3/T3;


AA3 = zeros(6);
AA3(6,3) = (1+aWIP3*hP3)/(T3*hA3);


%% define input

BD1=zeros(6,1);
BD1(2,1)=(ai1*hA1 - (1+aWIP1*hP1))/(T1*hA1);



BA1=zeros(6,1);
BA1(2,1)=(1+aWIP1*hP1)/(T1*hA1);


%% right hand side
rhs = A0*xx(:,1)...
      + AD1*xx(:, 2) + AP1*xx(:, 3) + AT1*xx(:, 4) + AA1*xx(:, 5)...
      + BD1*uu(:, 2) +                               BA1*uu(:, 5)...
      + AD2*xx(:, 6) + AP2*xx(:, 7) + AT2*xx(:, 8) + AA2*xx(:, 9)...
      + AD3*xx(:,10) + AP3*xx(:,11) + AT3*xx(:,12) + AA3*xx(:,13);


end

