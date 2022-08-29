% One big element
n = 32;
[x, y] = chebpts2(n, n, [-1 1 -0.5 0.5]);
z = 1+0*x;
nrm = sqrt(x.^2 + y.^2 + z.^2);
x = x ./ nrm;
y = y ./ nrm;
z = z ./ nrm;
dom1 = surfacemesh({x}, {y}, {z});
 

% Two glued elements
n = 16;
[xL, yL] = chebpts2(n, n, [-1 0 -0.5 0.5]);
[xR, yR] = chebpts2(n, n, [ 0 1 -0.5 0.5]);
zL = 1+0*xL;
zR = 1+0*xR;
nrmL = sqrt(xL.^2 + yL.^2 + zL.^2);
nrmR = sqrt(xR.^2 + yR.^2 + zR.^2);
xL = xL ./ nrmL;  xR = xR ./ nrmR;
yL = yL ./ nrmL;  yR = yR ./ nrmR;
zL = zL ./ nrmL;  zR = zR ./ nrmR;
dom2 = surfacemesh({xL;xR}, {yL;yR}, {zL;zR});

%% 
pdo = [];
pdo.lap = 1;
rhs1 = -2*spherefun.sphharm(1,0);
rhs1 = surfacefun(@(x,y,z) rhs1(x,y,z), dom1);
L1 = surfaceop(dom1,pdo,rhs1);
rhs2 = -2*spherefun.sphharm(1,0);
rhs2 = surfacefun(@(x,y,z) rhs2(x,y,z), dom2);
L2 = surfaceop(dom2,pdo,rhs2);

load('test.mat')
R1 = L1.patches{1}.R;
R2 = R;
 
leftIdx1  = 1:n-2;
downIdx1  = 1*(n-2)+1:2*(n-2);
upIdx1    = 2*(n-2)+1:3*(n-2);
downIdx2  = 3*(n-2)+1:4*(n-2);
upIdx2    = 4*(n-2)+1:5*(n-2);
rightIdx2 = 5*(n-2)+1:6*(n-2);

P = zeros(6*(n-2));
P(1:n-2, leftIdx1)              = eye(n-2);
P(1*(n-2)+1:2*(n-2), rightIdx2) = eye(n-2);
P(2*(n-2)+1:3*(n-2),downIdx1)   = eye(n-2);
P(3*(n-2)+1:4*(n-2),downIdx2)   = eye(n-2);
P(4*(n-2)+1:5*(n-2),upIdx1)     = eye(n-2);
P(5*(n-2)+1:6*(n-2),upIdx2)     = eye(n-2);
R2 = P * R2 * P';

[x0, ~, v0] = chebpts(n-2,[-1 0],1);
[x1, ~, v1] = chebpts(n-2,[0 1],1);
xh = [x0;x1];
vh = [v0;v1];
[x2, ~, v2] = chebpts(n-2,[-0.5, 0.5],1);
[x3, ~, v3] = chebpts(2*n-2,[-0.5, 0.5],1);
[x4, ~, v4] = chebpts(2*n-2,1);


C1 = barymat(x4,xh,vh);
C2 = barymat(x3,x2,v2);
CC = blkdiag(C2,C2,C1,C1);

C1 = barymat(xh,x4,v4);
C2 = barymat(x2,x3,v3);
CC1 = blkdiag(C2,C2,C1,C1);
%% 
norm(CC*R2*CC1 - R1)
