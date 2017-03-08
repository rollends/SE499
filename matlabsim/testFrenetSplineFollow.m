%% Select a few control points and generate appropriate polynomial to follow
clear all;
close all;
clc;

%Yc = [1, 2, 3, 4, 5, 4, 3, 2, 1; 1, 5, 4, 0, -2, -4, -3, -2, -1];
%Yc = [1, 2, 3, 4, 5; 1, 5, 4, 0, -2];
Yc = [1, 2, 3; 1, 5, 4];
Tc = ((1:size(Yc,2)) - 1) ./ (size(Yc,2) - 1);
pp = restrictedpolyfit(5, 1e-1, Yc);

controller = polysfcontrol(pp);


% Simulate Control System

% Initial Condition
X0 = [1, 1, pi/(2.5)];

global eta;
eta = -1;

options = odeset( 	'RelTol',1e-12, ...
                    'AbsTol', 1e-12, ...
                    'OutputFcn', @plotdriver, ...
                    ...'InitialStep', 1e-3, ...
                    ...'MaxStep', 1e-3, ...
                    'NormControl', 'on' );
odesolve = @ode45;
t_span = linspace(0, 20, (20-0)*50); % amount of time we will simulate for

% Translational velocity
v = 1;

figure;
curve = vpolyval(pp, linspace(0, 1, 100));
plot(curve(1,:), curve(2,:));
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(Yc(1,:), Yc(2,:), 'r*');
title('Output of differential drive robot');
grid on;

[t, x] = odesolve(@diffdrivesys, t_span, X0, options, v, controller);
disp('complete');

%% Save as video
driver_plot_video = VideoWriter('follow_SF_polynomial_deviate.mp4', 'MPEG-4');
driver_plot_video.FrameRate = 50;
open(driver_plot_video);

figure('position', [0 0 1024 1024]);
curve = vpolyval(pp, linspace(0, 1, 100));
plot(curve(1,:), curve(2,:));
xlabel('x_1(t)');
ylabel('x_2(t)');
xlim([1, 5]);
ylim([-5, 6]);
hold on;
plot(Yc(1,:), Yc(2,:), 'r*');
title('Output of Differential Drive Robot');
grid on;
for it = 1:length(t)
    plot(x(it,1), x(it,2), '.b');
    drawnow;
    writeVideo(driver_plot_video, getframe())
end
close(driver_plot_video);
clear('driver_plot_video');
disp('finished writing');