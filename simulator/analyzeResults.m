%% Load file (choose file name)
T = readtable('results_rrt_5.csv');

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
for i = 0:6
    x0 = 7 + 12 * i
    y0 = 5
    w = 5
    h = 20;
    fill([x0, x0+w, x0+w, x0],[y0, y0, y0+h, y0+h], 'black');
end
fill([0, 80, 80, 0],[50, 50, 60, 60], 'black');
fill([90, 100, 100, 90],[50, 50, 60, 60], 'black');
title('Robot World'); hold off;

%%
figure;
colormap lines;
for p = 1:rrtCount
    for i = 0:splineCount
        selSpline = (T.spline_ind == i) & (T.plan_ind == p);
        plot(T.xi_1(selSpline), T.xi_1(selSpline)); hold on;
    end
end
xlabel('Xi_1');
ylabel('Xi_2');
title('Linearized State'); hold off;