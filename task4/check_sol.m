clear all;
l = 5;
n = 20;
dom = surfacemesh.sphere(n, 2);

% build the true solution, which is spherical harmonics 
for i = 1:2*l+1
    
    sol = spherefun.sphharm(l,i-l-1);
    sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
    A(:,i) = sol;
    
end

% construct the operator 
pdo = [];
pdo.lap = 1;
pdo.b = l*(l+1);

tic

% compute the orthonormal eigenfunctions 
L = surfaceop(dom, pdo, 0);
Z = null(L);
Z = orth(orth(Z));

% compare coefficients
for i = 1: 2*l+1
    c = A\Z(:,i);
    Z1 = A*c;
    d = Z1 - Z(:,i);
    e(i) = norm(d);
end
toc