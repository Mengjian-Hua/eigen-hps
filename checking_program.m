% One big element
n = 10;
[x, y] = chebpts2(n, n, [-1 1 -0.5 0.5]);
z = 1+0*x;
dom1 = surfacemesh({x}, {y}, {z});
 

% Two glued elements
n = 10;
[xL, yL] = chebpts2(n, n, [-1 0 -0.5 0.5]);
[xR, yR] = chebpts2(n, n, [ 0 1 -0.5 0.5]);
zL = 1+0*xL;
zR = 1+0*xR;
dom2 = surfacemesh({xL;xR}, {yL;yR}, {zL;zR});

%% 

pdo = [];
pdo.lap = 1;

u_true_1 = surfacefun(@(x,y,z) x, dom1);
u_true_2 = surfacefun(@(x,y,z) x, dom2);

L1 = surfaceop(dom1,pdo,0);
L2 = surfaceop(dom2,pdo,0);

build(L1)
build(L2)


bc = @(x,y,z) x;
u1 = solve(L1,bc);
u2 = solve(L2,bc);

% bc1 = eta*L1.patches{1}.xyz(:,1);
% idx = abs(L1.patches{1}.xyz(:,1) -1)<1e-13;% idx([end-32 end]) = false;
% bc1(idx) = bc1(idx) + 1;
% idx = abs(L1.patches{1}.xyz(:,1) +1)<1e-13;% idx([end-2*32+1 end-32+1]) = false;
% bc1(idx) = bc1(idx) - 1;
% [xn1, ~, vn1] = chebpts(32, 2);
% [xn2, ~, vn2] = chebpts(32, 1);
% C  = barymat(xn2, xn1, vn1);
% CC = blkdiag(C, C, C, C);
% bc1 = CC*bc1;

% bc2 = eta*L2.patches{1}.xyz(:,1);
% idx = abs(L2.patches{1}.xyz(:,1) -1)<1e-13;% idx([end-32 end]) = false;
% bc2(idx) = bc2(idx) + 1;
% idx = abs(L2.patches{1}.xyz(:,1) +1)<1e-13;% idx([end-2*32+1 end-32+1]) = false;
% bc2(idx) = bc2(idx) - 1;
% bc2 = CC*bc2;

