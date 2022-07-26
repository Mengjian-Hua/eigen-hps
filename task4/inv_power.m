%% inverse power iteration for spherical harmonics
clear all;
l = 5; % number of eigenvalues we want to compute
n = 30; % number of mesh points
dom = surfacemesh.blob(n);

pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, 0);
tol = 1e-3; % tolerance used to stop the iteration

% build a random orthonormal basis 
for i = 1:l
    Q(:,i) = surfacefun(@(x,y,z) sin(i*x)+cos(i*y), dom);
end
Q = orth(orth(Q));

% first iteration
tic
for i = 1:l
    L = updateRHS(L, Q(:,i));
    Z(:,i) = L.solve();
end

[Q,R] = qr(Z);
num_iter = 1; % number of iterations
R_previous = zeros(size(R));

% subsequent iterations
while(trace(abs(R-R_previous))>tol)
    
    for i = 1:l
        L = updateRHS(L, Q(:,i));
        Z(:,i) = L.solve();
    end
    
    R_previous = R;
    [Q,R] = qr(Z);
    num_iter = num_iter+1;

    
end
toc
