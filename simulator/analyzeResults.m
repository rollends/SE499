%% Load file (choose file name)
T = readtable('results.csv');

splineCount = max(T.spline_ind);

%%
figure;
colormap lines;
plot(T.x_1, T.x_2); hold on;
for i = 0:splineCount
    selSpline = (T.spline_ind == i);
    plot(T.sigma_1(selSpline), T.sigma_2(selSpline));
end
xlabel('World X');
ylabel('World Y');
title('Robot World'); hold off;

figure;
colormap lines;
for i = 0:splineCount
    selSpline = (T.spline_ind == i);
    plot(T.xi_1(selSpline), T.xi_1(selSpline)); hold on;
end
xlabel('Xi_1');
ylabel('Xi_2');
title('Linearized State'); hold off;