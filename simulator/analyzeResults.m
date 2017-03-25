%% Load file (choose file name)
T = readtable('results_rrt_sfcontrol.csv');

splineCount = max(T.spline_ind);
rrtCount = max(T.plan_ind);

figure;
colormap lines;
figBot = plot(T.x_1, T.x_2, '--'); 
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
legend([figBot, figTransfer], 'Robot Path', 'RRT (planning) Events');
for i = 0:5
    x0 = 10 + 10 * i;
    y0 = 5;
    w = 5;
    h = 20;
    fill([x0, x0+w, x0+w, x0],[y0, y0, y0+h, y0+h], 'black');
end
fill([0, 80, 80, 0],[50, 50, 60, 60], 'black');
fill([90, 100, 100, 90],[50, 50, 60, 60], 'black');
title('Robot World'); hold off;

figure(2);
plot(T.time, T.xi_1); 
ylabel('Error in xi_1');
hold on;
%plot(T.time, sqrt( (T.sigma_1-T.x_1) .^ 2 + (T.sigma_2-T.x_2) .^ 2 - (T.xi_1 .^ 2) ));
yyaxis right;
plot(T.time, T.delta_xi_2);
title('Orthogonal Tracking Error');
ylabel('Error in xi_2');
xlabel('Time (s)');
hold off;

figure(3);
for p = 1:rrtCount
    transferPoints = (T.plan_ind == p);
    plot(T.xi_1(transferPoints), T.xi_2(transferPoints), '--'); hold on;
end
D = 0.1^4;
v = 1;
b = D/(v^4);
a = 2*sqrt(D) / (v^2);
k = v^2 / D^0.25;
V = @(x,y) x.^4 + a*(x.^2).*(y.^2) + b * y.^4 - D;
fimplicit(V);
title('Linearized State');
xlabel('xi_1');
ylabel('xi_2');
hold off;

%% Campus
T = readtable('results_rrt_sfcontrol_campus.csv');

splineCount = max(T.spline_ind);
rrtCount = max(T.plan_ind);

Map = imread('../report/campus_map.PNG');
Map = Map(end:-1:1,:,:);

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
plot(T.time, T.xi_1); 
ylabel('Error in xi_1');
hold on;
yyaxis right;
plot(T.time, T.delta_xi_2);
title('Orthogonal Tracking Error');
ylabel('Error in xi_2');
xlabel('Time (s)');
hold off;

figure(3);
for p = 1:rrtCount
    transferPoints = (T.plan_ind == p);
    plot(T.xi_1(transferPoints), T.xi_2(transferPoints), '--'); hold on;
end
D = 0.1^4;
v = 1;
b = D/(v^4);
a = 2*sqrt(D) / (v^2);
k = v^2 / D^0.25;
V = @(x,y) x.^4 + a*(x.^2).*(y.^2) + b * y.^4 - D;
fimplicit(V);
title('Linearized State');
xlabel('xi_1');
ylabel('xi_2');
hold off;
