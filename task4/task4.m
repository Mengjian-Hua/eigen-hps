%% Plot the monitor function for the Laplace-Beltrami operator

n = 30;
dom = surfacemesh.sphere(n, 2);

%%
lambda_list = -0.1:-0.2:-21;
tic
parfor i = 1:length(lambda_list)
    lambda = lambda_list(i);
    pdo = [];
    pdo.lap = 1;
    monitor(i) = eval_monitor(lambda,pdo,dom);
end
toc
plot(lambda_list,monitor,"linewidth",2)