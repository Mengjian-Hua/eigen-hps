%% inverse power iteration for spherical harmonics
clear all;
l = 16;
n = 20;
dom = surfacemesh.sphere(n, 2);

pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, 0);
tol = 1e-3;

for i = 1:l
    Q(:,i) = surfacefun(@(x,y,z) sin(i*x)+cos(i*y), dom);
end
Q = orth(orth(Q));

tic
for i = 1:l
    L = updateRHS(L, Q(:,i));
    Z(:,i) = L.solve();
end

[Q,R] = qr(Z);
num_iter = 1;
R_previous = zeros(size(R));

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
