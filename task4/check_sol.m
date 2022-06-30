clear all;
l = 5;
n = 20;
dom = surfacemesh.sphere(n, 2);

for i = 1:2*l+1
    
    sol = spherefun.sphharm(l,i-l-1);
    sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
    A(:,i) = sol;
    
end

pdo = [];
pdo.lap = 1;
pdo.b = l*(l+1);
tic
L = surfaceop(dom, pdo, 0);
Z = null(L);
Z = orth(Z);
for i = 1: 2*l+1
    c = A\Z(:,i);
    Z1 = A*c;
    d = Z1 - Z(:,i);
    e(i) = norm(d);
end
toc