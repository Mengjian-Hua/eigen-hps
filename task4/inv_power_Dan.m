%% inverse power iteration for spherical harmonics
clear Q R Z
rng(0)
l = 4; % Number of eigenvalues to compute
n = 16;
bc = 0;
%dom = surfacemesh.cyclide(n, 16, 16);
%dom = surfacemesh.blob(n, 2);
%dom = surfacemesh.sphere(n, 2);
%dom = surfacemesh.ellipsoid(n, [4.5 6 3], 2);
dom = surfacemesh.mobius(n, 16, 4);
%dom = surfacemesh.cube(n, 2);
%dom = resample(surfacemesh.fromRhino('models/cow.csv', 8), n);

pdo = [];
pdo.lap = 1;
pdo.b = 300; % Shift
L = surfaceop(dom, pdo);
L.rankdef = true;
tol = 1e-8;

% Build a random orthonormal basis
for i = 1:l
    Q(:,i) = surfacefun(@(x,y,z) rand(size(x))-1, dom);
end
Q = orth(Q);

tic
L.rhs = Q;
Z = L.solve(bc);
[Q, R] = qr(Z);
num_iter = 1;
R_previous = zeros(size(R));
while ( num_iter == 1 || max(abs(R-R_previous), [], 'all') > tol )
    num_iter
    L.rhs = Q;
    Z = L.solve(bc);
    R_previous = R;
    [Q, R] = qr(Z);
    num_iter = num_iter+1;
end
r = diag(R);
toc
