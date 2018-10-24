function [ tauVec ] = supplyChainDelays( xx,alpha,~ )
% 
%
approx=300^2;

TransportDelay1=(sqrt(approx-approx*exp(-xx(1,1)/approx)) + 1);
TransportDelay2=(sqrt(approx-approx*exp(-xx(3,1)/approx)) + 1);
TransportDelay3=(sqrt(approx-approx*exp(-xx(5,1)/approx)) + 1);

tauD1 = alpha(1);
tauP1 = alpha(1) + alpha(2);
tauT1 = alpha(1) + alpha(2) + TransportDelay1;
tauA1 = alpha(1) + 11;

tauD2 = alpha(3);
tauP2 = alpha(3) + alpha(4);
tauT2 = alpha(3) + alpha(4) + TransportDelay2;
tauA2 = alpha(3) + 12;

tauD3 = alpha(5);
tauP3 = alpha(5) + alpha(6);
tauT3 = alpha(5) + alpha(6) + TransportDelay3;
tauA3 = alpha(5) + 13;

tauVec=[tauD1, tauP1, tauT1, tauA1,...
        tauD2, tauP2, tauT2, tauA2,...
        tauD3, tauP3, tauT3, tauA3]';

end

 
