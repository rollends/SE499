function dx = diffdrivesys( t, x, v, controller )
%DIFFDRIVESYS Non-linear system modelling a basic differential drive
%system.
%   Standard model. Use odeXX family to solve. u is control input.

x1 = x(1);
x2 = x(2);
x3 = x(3);

% controller
u = controller(x);

% the system

f = [v*cos(x3);v*sin(x3);0]; g = [0;0;1];
dx =  f + g*u;

