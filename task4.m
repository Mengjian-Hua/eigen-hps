% Laplace-Beltrami

n = 20;
dom = surfacemesh.sphere(n, 2);

%%
lambda_list = -0.1:-0.2:-20;
tic
for i = 1:length(lambda_list)
    lambda = lambda_list(i);
    pdo = [];
    pdo.lap = 1;
    monitor(i) = eval_monitor(lambda,pdo,dom);
end
toc
plot(lambda_list,monitor,"linewidth",2)