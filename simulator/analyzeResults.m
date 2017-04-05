%% Campus
T = readtable('results_rrt_sylvester_campus.csv');

splineCount = max(T.spline_ind);
rrtCount = max(T.plan_ind);

Map = imread('../report/images/campus_map.PNG');
Map = Map(end:-1:1,:,:) * 1.75;

figure;
colormap lines;
image('CData', Map); hold on;
figBot = plot(T.x_1, T.x_2, '.', 'MarkerSize', 5); 
hold on;
for p = 1:rrtCount
    for i = 0:splineCount
        selSpline = (T.spline_ind == i) & (T.plan_ind == p);
        plot(T.sigma_1(selSpline), T.sigma_2(selSpline));
    end
    
    transferPoints = (diff(T.plan_ind) ~= 0);
    figTransfer = scatter(T.x_1(transferPoints), T.x_2(transferPoints), '*');
end
xlabel('World X');
ylabel('World Y');
xlim([0, size(Map,2)]);
ylim([0, size(Map,1)]);
legend([figBot, figTransfer], 'Robot Path', 'RRT (planning) Events');
title('Robot World'); hold off;

figure(2);
plot(T.time, -sqrt(10) * T.xi_1 - T.xi_2 * 2 * (10^0.25)); 
ylabel('u');
title('Control Effort');
xlabel('Time (s)');
hold off;


%% Campus for Tracking Controller
T = readtable('results_rrt_trackcontrol_campus_mod.csv');

splineCount = max(T.spline_ind);
rrtCount = max(T.plan_ind);

Map = imread('../report/images/campus_map.PNG');
Map = Map(end:-1:1,:,:) * 1.75;

figure;
colormap lines;
image('CData', Map); hold on;
figBot = plot(T.x_1, T.x_2, '.', 'MarkerSize', 5); 
hold on;
for p = 1:rrtCount
    for i = 0:splineCount
        selSpline = (T.spline_ind == i) & (T.plan_ind == p);
        plot(T.sigma_1(selSpline), T.sigma_2(selSpline));
    end
    
    transferPoints = (diff(T.plan_ind) ~= 0);
    figTransfer = scatter(T.x_1(transferPoints), T.x_2(transferPoints), '*');
end
xlabel('World X');
ylabel('World Y');
xlim([0, size(Map,2)]);
ylim([0, size(Map,1)]);
legend([figBot, figTransfer], 'Robot Path', 'RRT (planning) Events');
title('Robot World'); hold off;

figure(2);
plot(T.time, sqrt((T.x_1 + 1*cos(T.x_3) - T.sigma_1).^2 + (T.x_2 + 1*sin(T.x_3) - T.sigma_2).^2))
title('Tracking Error');
ylabel('Error');
xlabel('Time (s)');
hold off;
