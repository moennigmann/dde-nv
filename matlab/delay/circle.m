function [ x,y ] = circle( radius, x, y )
% gives back the x-y-coordinates of a circle with defined radius and center
% at x,y

phi=0:0.01:2*pi;

x=radius*sin(phi)+x;
y=radius*cos(phi)+y;

end

