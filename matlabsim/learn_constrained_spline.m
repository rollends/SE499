%% Select a few control points
clear all;
close all;
clc;
digits(50);
useVPA = false;
if useVPA 
    type = @vpa;
else 
    type = @double;
end

Yc = [1, 2, 3, 4, 5; 1, 5, 4, 0, -2];
Yc = [1, 2, 3, 4, 5, 4, 3, 2, 1; 1, 5, 4, 0, -2, -4, -3, -2, -1];
Tc = ((1:size(Yc,2)) - 1) ./ (size(Yc,2) - 1);

%% Generate a BUNCH of points in between these control points
Y = Yc(:,1);
N = 12;
Yc = type(Yc);
SampleRate = 20;
for i = 1:(size(Yc, 2)-1)
    samples = vertcat( linspace(Yc(1,i),Yc(1,i+1),SampleRate),...
                       linspace(Yc(2,i),Yc(2,i+1),SampleRate) );
    Y = horzcat(Y, samples(:,2:end));
end

%% Construct Least Squares Optimization Problem
% Would like to find coefficients for 30th order polynomial (31x2 coeffs).
C = type(zeros(numel(Y), (N+1) * 2));
d = type(zeros(numel(Y), 1));

index = 1;
% Go through every coordinate
for xi = 1:size(Y, 1)
    % Go through every point
    for pi = 1:size(Y, 2)
        base = (xi-1)*(N+1);
        lambda = type( (pi - 1) / (size(Y,2) - 1) );
        polyLambda = lambda .^ (N:-1:0);
        
        C(index, (base+1):(base+(N+1))) = polyLambda;
        d(index) = Y(xi, pi);

        % Next row of C
        index = index + 1;
    end    
end
rank(C)

%% Perform basic least squares
p1 = pinv(C)*d;
T = linspace(0, 1, size(Y, 2));
% Ys = [polyval(p1(1:(N+1)), T); polyval(p1(N+2:end), T)];
% Yerror = (Ys - Y).^2;
% Yerror = sqrt( sum(Yerror, 1) );
% Yerror = sum(Yerror);

%% Construct bounds for control points
% Want < 1e-1 error for control points. (here constrained by L1 norm)
A = type(zeros(2 * numel(Yc), (N+1) * 2));
b = type(zeros(2 * numel(Yc), 1));
delta = 0.2;

index = 1;
% Go through every coordinate
for xi = 1:size(Yc, 1)
    % Go through every point
    for pi = 1:size(Yc, 2)
        base = (xi-1)*(N+1);
        lambda = type( (pi - 1) / (size(Yc,2) - 1) );
        polyLambda = lambda .^ (N:-1:0);
 
        % Enforce upper bound
        A(index, (base+1):(base+(N+1))) = polyLambda;
        b(index) = Yc(xi, pi) + delta;

        % Next row of C
        index = index + 1;
        
        % Enforce lower bound
        A(index, (base+1):(base+(N+1))) = -polyLambda;
        b(index) = -Yc(xi, pi) + delta;
        
        % Next row of C
        index = index + 1;
    end
end

%% Enforce End Point Conditions (Fixed End Points)
Aeq = type(zeros(4, (N+1) * 2));
beq = type(zeros(4, 1));

Aeq(1, 1:(N+1)) = circshift(eye(1, N+1), -1);
Aeq(2, (N+2):end) = circshift(eye(1, N+1), -1);
beq(1:2) = Yc(:,1);

Aeq(3, 1:(N+1)) = ones(1, N+1);
Aeq(4, (N+2):end) = ones(1, N+1);
beq(3:4) = Yc(:,end);

%% Perform Constrained Least Squares
option = optimoptions(  @lsqlin, ...
                        'Algorithm', 'interior-point', ...
                        'MaxIterations', 1000, ...
                        'ConstraintTolerance', delta * 1e-1);
[p2,error] = lsqlin(double(C), double(d), ... Least Squares System
            double(A), double(b), ... Linear Constraints (1-Norm Ball around Control Points)
            double(Aeq), double(beq), ... Equality Constraints (ignored)
            [], ...-1e3 * ones(size(p1)), ... Lower Bound
            [], ...1e3 * ones(size(p1)), ... Upper Bound
            double(p1), ...
            option);
        
%% Perform Constrained Least Squares by QuadProg
% Cs = C;
% Q = (Cs')*Cs;
% q = -(Cs')*d;
% 
% option = optimoptions(  @quadprog, ...
%                         'Algorithm', 'interior-point-convex' );
% p2 = quadprog(  double(Q), double(q), ... Least Squares Matrices
%                 double(A), double(b), ... Linear Constraints
%                 double(Aeq), double(beq), ... Equality Constraints
%                 [], ... Lower Bound
%                 [], ... Upper Bound
%                 [], ... Initializer
%                 option );

%% Plot both
T = linspace(0, 1, 1000);
p1 = double(p1);
p2 = double(p2);
C1 = vertcat(polyval(p1(1:(N+1)), T), polyval(p1((N+1+1):end), T));
C2 = vertcat(polyval(p2(1:(N+1)), T), polyval(p2((N+1+1):end), T));

close all;
plot(Yc(1,:), Yc(2,:), 'bo'); hold on;
plot(C1(1,:), C1(2,:), 'r-');
plot(polyval(p1(1:(N+1)), Tc), ...
     polyval(p1((N+1+1):end), Tc),...
     'r*');
plot(C2(1,:), C2(2,:), 'g-');
plot(polyval(p2(1:(N+1)), Tc), ...
     polyval(p2((N+1+1):end), Tc),...
     'g*');
hold off;
%plot(Y(1,:), Y(2,:), 'b.');
