% One big element
n = 32;
[x, y] = chebpts2(n, n, [-1 1 -0.5 0.5]);
z = 1+0*x;
% nrm = sqrt(x.^2 + y.^2 + z.^2);
% x = x ./ nrm;
% y = y ./ nrm;
% z = z ./ nrm;
dom1 = surfacemesh({x}, {y}, {z});
 

% Two glued elements
n = 32;
[xL, yL] = chebpts2(n, n, [-1 0 -0.5 0.5]);
[xR, yR] = chebpts2(n, n, [ 0 1 -0.5 0.5]);
zL = 1+0*xL;
zR = 1+0*xR;
% nrmL = sqrt(xL.^2 + yL.^2 + zL.^2);
% nrmR = sqrt(xR.^2 + yR.^2 + zR.^2);
% xL = xL ./ nrmL;  xR = xR ./ nrmR;
% yL = yL ./ nrmL;  yR = yR ./ nrmR;
% zL = zL ./ nrmL;  zR = zR ./ nrmR;
dom2 = surfacemesh({xL;xR}, {yL;yR}, {zL;zR});

%% 

pdo = [];
pdo.lap = 1;

u_true_1 = surfacefun(@(x,y,z) x, dom1);
u_true_2 = surfacefun(@(x,y,z) x, dom2);
eta = 2i;

% ii = abs(x) < 1 & abs(y) < 0.5;     % Interior indices.
% ee = ~ii;

% boundary_1 = eta*x;
% boundary_1(ii) = 0;
% boundary_1 = nonzeros(boundary_1(:));
% boundary_1(1:32) = boundary_1(1:32) - 1;
% boundary_1(93:end) = boundary_1(93:end) + 1;
% 
% x2 = [xL xR];
% y2 = [yL yR];
% boundary_2 = eta*x2;
% ii = abs(x2) < 1 & abs(y2) < 0.5;     % Interior indices.
% ee = ~ii;
% boundary_2 = boundary_2(ee);
% boundary_2 = boundary_2(:);
% boundary_2(1:16) = boundary_2(1:16) - 1;
% boundary_2(77:92) = boundary_2(77:92) + 1;


L1 = surfaceop(dom1,pdo,0);
L2 = surfaceop(dom2,pdo,0);

build(L1)
build(L2)

bc1 = eta*L1.patches{1}.xyz(:,1);
idx = abs(L1.patches{1}.xyz(:,1) -1)<1e-13;% idx([end-32 end]) = false;
bc1(idx) = bc1(idx) + 1;
idx = abs(L1.patches{1}.xyz(:,1) +1)<1e-13;% idx([end-2*32+1 end-32+1]) = false;
bc1(idx) = bc1(idx) - 1;
% [xn1, ~, vn1] = chebpts(32, 2);
% [xn2, ~, vn2] = chebpts(32, 1);
% C  = barymat(xn2, xn1, vn1);
% CC = blkdiag(C, C, C, C);
% bc1 = CC*bc1;

bc2 = eta*L2.patches{1}.xyz(:,1);
idx = abs(L2.patches{1}.xyz(:,1) -1)<1e-13;% idx([end-32 end]) = false;
bc2(idx) = bc2(idx) + 1;
idx = abs(L2.patches{1}.xyz(:,1) +1)<1e-13;% idx([end-2*32+1 end-32+1]) = false;
bc2(idx) = bc2(idx) - 1;
% bc2 = CC*bc2;

u1 = solve(L1,bc1);
u2 = solve(L2,bc2);

