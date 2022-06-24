%% Finding eigenvalues of the Laplace-Beltrami operator (example)
n = 10;
dom = surfacemesh.sphere(n, 2);
pdo = [];
pdo.lap = 1;
% initial guess
[V,D] = eig_finder(dom,pdo);