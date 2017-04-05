%% Create Curve
clear all;
close all;
clc;

Yc = [1, 2, 3, 4, 5; 1, 5, 4, -3, -3];
Tc = ((1:size(Yc,2)) - 1) ./ (size(Yc,2) - 1);
pp = restrictedpolyfit(5, 1, Yc);

% Initial Condition
X0 = [1, 1, 0];

% Solver
options = odeset( 	'RelTol',1e-7, ...
                    'AbsTol', 1e-7, ...
                    'OutputFcn', @plotdriver );
odesolve = @ode45;
t_span = [0, 15];%linspace(0, 15, (15-0)*500);

%% Sylvester Follower
close all;

figure(1);
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

figure(2);
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
% close all;
% 
% figure(3);
% curve = vpolyval(pp, linspace(0, 1, 500));
% plot(curve(1,:), curve(2,:));
% xlabel('x_1(t)');
% ylabel('x_2(t)');
% xlim([1, 5]);
% ylim([-5, 6]);
% hold on;
% plot(Yc(1,:), Yc(2,:), 'r*');
% title('Output of differential drive robot');
% grid on;
% 
% robot3 = TrackingRobot(pp);
% [tF, xF] = odesolve(@(t, x) robot3.stepRobot(t, x), t_span, X0, options);
% 
% disp('Completed Tracking Follower');

%% Comparitive Plots

figure(4);
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
subplot(1, 2, 1);
plot(robot.solveXi(1, :), robot.solveXi(2,:), 'r-');
xlabel('\xi_1(t)');
ylabel('\xi_2(t)');
title('Sylvester');

subplot(1, 2, 2);
plot(robot2.solveXi(1, :), robot2.solveXi(2,:), 'b-');
xlabel('\xi_1(t)');
ylabel('\xi_2(t)');
title('Serret-Frenet');

figure;
plot(robot.solveT,sqrt(sum((robot.solveX(1:2,:) - robot.solveS) .^ 2, 1))); hold on;
plot(robot2.solveT,sqrt(sum((robot2.solveX(1:2,:) - robot2.solveS) .^ 2, 1)));

figure;
plot(robot.solveT,robot.solveU); hold on;
plot(robot2.solveT,robot2.solveU);

%%
% 
% figure;
% curve = vpolyval(pp, linspace(0, 1, 500));
% plot(curve(1,:), curve(2,:), '-');
% xlabel('x_1(t)');
% ylabel('x_2(t)');
% xlim([1, 5]);
% ylim([-5, 6]);
% hold on;
% plot(robot3.solveX(1, :) + 0.05*cos(robot3.solveX(3, :)), robot3.solveX(2,:) + 0.05*sin(robot3.solveX(3, :)), '--');
% xlabel('x_1(t)');
% ylabel('x_2(t)');
% title('Tracking of Tracking Controlled Robot');
