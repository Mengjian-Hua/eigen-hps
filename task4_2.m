%% Finding the root of the monitor function

% Laplace-Beltrami

n = 20;
dom = surfacemesh.sphere(n, 2);

% initial guess
lambda = -1;
diff = 1e-2;
lambda2 = lambda + diff;

f1 = eval_monitor(lambda,pdo,dom);

num_eval = 1;
tol = 1e-4;

while (f1>tol)
    f2 = eval_monitor(lambda2,pdo,dom);
    df = (f2 - f1)/diff;
    lambda = lambda - f1/df;
    lambda2 = lambda + diff;
    f1 = eval_monitor(lambda,pdo,dom);
    num_eval = num_eval + 2;
end