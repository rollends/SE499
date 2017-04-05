% Chosen bounds
M = 20;
v = 1;

% Parameters of the lyapunov function
D = M^4;
b = D/(v^4);
a = 2*sqrt(D) / (v^2);

% The gain required to ensure stability (only g(1) matters, I've chosen g(2) here
% arbitrarily)
g = [-(v/M)^2, -v / M * sqrt(2)];

% The Lyapunov function, offset by D so that we can also use it to plot the open,
% connected set Omega
V = @(x,y) x.^4 + a*(x.^2).*(y.^2) + b * y.^4 - D;
V = @(x,y) (x.^2 + (M/v)^2 .* y.^2).^2 - D; % this definition is just a re-arranged version of above.

% Form a grid to work with and calculate the flow on the grid.
[xi1,xi2] = meshgrid(-(D^(0.25)):0.5:(D^0.25), -1:0.05:1);
dxi1 = xi2;
dxi2 = g(1) * xi1 + g(2)*xi2;

% Plot the lyapunov function over the set Omega and with the flow
figure(1);
surf(xi1, xi2, V(xi1, xi2) + D, 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on;
quiver(xi1, xi2, dxi1, dxi2);
fimplicit(V);
xlabel('\xi_1');
ylabel('\xi_2');
zlabel('Energy');

% Plot just the flow and the open set
figure(2);
quiver(xi1, xi2, dxi1, dxi2); hold on;
fimplicit(V);
xlabel('\xi_1');
ylabel('\xi_2');