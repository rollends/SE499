%% Create Curve
clear all;
close all;
clc;

Yc = [1, 2, 3, 4, 5; 1, 5, 4, -3, -3];
Tc = ((1:size(Yc,2)) - 1) ./ (size(Yc,2) - 1);
pp = restrictedpolyfit(5, 1, Yc);

% Initial Condition
X0 = [1, 1, pi/4];

% Solver
options = odeset( 	'RelTol',1e-12, ...
                    'AbsTol', 1e-12, ...
                    'OutputFcn', @plotdriver );
odesolve = @ode45;
t_span = linspace(0, 15, (15-0)*200);

%% Sylvester Follower
close all;

figure;
curve = vpolyval(pp, linspace(0, 1, 500));
plot(curve(1,:), curve(2,:));
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(Yc(1,:), Yc(2,:), 'r*');
title('Output of differential drive robot');
grid on;

robot = SylvesterPathFollower(pp);

[tS, xS] = odesolve(@(t, x) robot.stepRobot(t, x), t_span, X0, options);
disp('Completed Sylvester Follower');

%% Frame Follower
close all;

figure;
curve = vpolyval(pp, linspace(0, 1, 500));
plot(curve(1,:), curve(2,:));
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(Yc(1,:), Yc(2,:), 'r*');
title('Output of differential drive robot');
grid on;

robot2 = FramePathFollower(pp);
[tF, xF] = odesolve(@(t, x) robot2.stepRobot(t, x), t_span, X0, options);

disp('Completed Frame Follower');

%% Tracking Follower
close all;

figure;
curve = vpolyval(pp, linspace(0, 1, 500));
plot(curve(1,:), curve(2,:));
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(Yc(1,:), Yc(2,:), 'r*');
title('Output of differential drive robot');
grid on;

robot3 = TrackingRobot(pp);
[tF, xF] = odesolve(@(t, x) robot3.stepRobot(t, x), t_span, X0, options);

disp('Completed Tracking Follower');

%% Comparitive Plots
figure;
curve = vpolyval(pp, linspace(0, 1, 500));
plot(curve(1,:), curve(2,:), 'b-');
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(robot.solveX(1,:), robot.solveX(2,:), 'r--');
plot(robot2.solveX(1,:), robot2.solveX(2,:), 'g--');
title('Track of Robots');

figure;
plot(robot.solveXi(1, :), robot.solveXi(2,:), 'r-');
xlabel('xi_1(t)');
ylabel('xi_2(t)');
title('Linearized State of Robot (Sylvester)');

figure;
plot(robot2.solveXi(1, :), robot2.solveXi(2,:), 'b-');
xlabel('xi_1(t)');
ylabel('xi_2(t)');
title('Linearized State of Robot (Frame Based)');

%%

figure;
curve = vpolyval(pp, linspace(0, 1, 500));
plot(curve(1,:), curve(2,:), 'b-');
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(robot3.solveX(1, :) + 0.05*cos(robot3.solveX(3, :)), robot3.solveX(2,:) + 0.05*sin(robot3.solveX(3, :)), 'b.');
xlabel('xi_1(t)');
ylabel('xi_2(t)');
title('Linearized State of Robot (Frame Based)');
