function [ F ] = Syncmanifold_3L_Setup1_fkt(var,pump,omega)
%   Berechnung von Allen Systemgrößen am Steadystate für 3 symmetrisch
%   gekoppelte Laserdioden   29.6.17   Jens Müller

x(1:9)=var(1:9);
xtau(1:9)=var(1:9);
alpha(1:3)=pump(1:3);
p(1)=omega;


F=[p(1) * x(2) + 0.5e0 * x(3) * x(1) - 0.20e1 * x(3) * x(2) + 0.5e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.5e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(2);
  -p(1) * x(1) + 0.20e1 * x(3) * x(1) + 0.5e0 * x(3) * x(2) - 0.5e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.5e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(2);
  0.5e-2 * alpha(1) - 0.5e-2 * x(3) - 0.5e-2 * (x(3) + 0.1e1) * (x(1)^2 + x(2)^2);
  p(1) * x(5) + 0.5e0 * x(6) * x(4) - 0.20e1 * x(6) * x(5) + 0.5e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.5e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(2);
  -p(1) * x(4) + 0.20e1 * x(6) * x(4) + 0.5e0 * x(6) * x(5) - 0.5e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.5e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(2);
  0.5e-2 * alpha(2) - 0.5e-2 * x(6) - 0.5e-2 * (x(6) + 0.1e1) * (x(4)^2 + x(5)^2);
  p(1) * x(8) + 0.5e0 * x(9) * x(7) - 0.20e1 * x(9) * x(8) + 0.5e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.5e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(2);
  -p(1) * x(7) + 0.20e1 * x(9) * x(7) + 0.5e0 * x(9) * x(8) - 0.5e-2 * sin(0.100e3 * p(1) + 0.2e1) * xtau(1) + 0.5e-2 * cos(0.100e3 * p(1) + 0.2e1) * xtau(2);
  0.5e-2 * alpha(3) - 0.5e-2 * x(9) - 0.5e-2 * (x(9) + 0.1e1) * (x(7)^2 + x(8)^2)];
end

