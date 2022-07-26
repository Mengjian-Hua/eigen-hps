%% Using rectangular ideas
rng(0)
u = randnfun2(1);
f = lap(u);

n = 10;
nbdy = n;
kindFrom = 1;
kindTo   = 1;
nu = 1;
eta = 2i;
numIntPts = (n-2)^2;
numBdyPts = 4*nbdy;

[xx, yy] = chebpts2(n);
bb = (xx == -1 | xx == 1 | yy == -1 | yy == 1);
ii = ~bb;
[x0, ~, v0] = chebpts(n,   kindFrom);
[x2, ~, v2] = chebpts(n-2, kindTo);
[xx0, yy0] = meshgrid(x0);
[xx2, yy2] = meshgrid(x2);

xn1 = chebpts(nbdy, 1);
xxb = [ -1+0*xn1 ; 1+0*xn1 ; xn1 ; xn1 ];
yyb = [ xn1 ; xn1 ; -1+0*xn1 ; 1+0*xn1 ];

B = zeros(numBdyPts, n^2);
for k = 1:numBdyPts
    B(k,:) = kron(barymat(xxb(k), x0, v0), barymat(yyb(k), x0, v0));
end

I = eye(n);
D = diffmat(n, 1, ['chebkind' num2str(kindFrom)]);
Dx = kron(D, I);
Dy = kron(I, D);

% Normal derivative on a square
N = zeros(numBdyPts, n^2);
N(1:nbdy,:)          = B(1:nbdy,:)          * -Dx;
N(nbdy+1:2*nbdy,:)   = B(nbdy+1:2*nbdy,:)   *  Dx;
N(2*nbdy+1:3*nbdy,:) = B(2*nbdy+1:3*nbdy,:) * -Dy;
N(3*nbdy+1:4*nbdy,:) = B(3*nbdy+1:4*nbdy,:) *  Dy;

S02 = barymat(x2, x0, v0);
D2 = diffmat([n-2 n], 2, ['chebkind' num2str(kindFrom)], ['chebkind' num2str(kindTo)]);
L = kron(D2, S02) + kron(S02, D2);
ff = f(xx2, yy2);
uu = u(xx0, yy0);
gg = (nu*N + eta*B) * uu(:);

X = [L ; nu*N + eta*B] \ [zeros(numIntPts, numBdyPts) ff(:); eye(numBdyPts) zeros(numBdyPts, 1)];
Y = (nu*N - eta*B) * X(:,1:numBdyPts);

[xn1, ~, vn1] = chebpts(nbdy, 1);
[xn2, ~, vn2] = chebpts(n-2, 1);
C  = barymat(xn2, xn1, vn1);
C1 = barymat(xn1, xn2, vn2);
CC = blkdiag(C, C, C, C);
CC1 = blkdiag(C1, C1, C1, C1);
Y = CC * Y * CC1;

norm(uu(:) - X*[gg(:) ; 1])
cond(Y)
figure(1), plot(eig(Y), 'o')
shg

%% Removing the corners
rng(0)
u = randnfun2(1);
f = lap(u);

n = 30;
eta = 2i;

[xx, yy] = chebpts2(n);
bb = (xx == -1 | xx == 1 | yy == -1 | yy == 1);
bb([1 end], [1 end]) = false;
ii = ~bb;
ii([1 end], [1 end]) = false;
nc = ii | bb;
numBdyPts = 4*(n-2);
numIntPts = (n-2)^2;

I = eye(n);
D = diffmat(n, 1);
Dx = kron(D, I); Dx = Dx(bb,:);
Dy = kron(I, D); Dy = Dy(bb,:);
ibc = 3*(n-2)+1;
leftIdx  = 1:n-2;             leftIdxG  = 1:n-2;
rightIdx = n-1:2*(n-2);       rightIdxG = ibc:4*(n-2);
downIdx  = 2*(n-2)+1:3*(n-2); downIdxG  = n-1:2:ibc-2;
upIdx    = 3*(n-2)+1:4*(n-2); upIdxG    = n:2:ibc-1;

% Normal derivative on a square
N = zeros(4*(n-2), n^2);
N(leftIdx,:)  = -Dx(leftIdxG,:);
N(rightIdx,:) =  Dx(rightIdxG,:);
N(downIdx,:)  = -Dy(downIdxG,:);
N(upIdx,:)    =  Dy(upIdxG,:);

I = eye(n);
II = kron(I, I);
D2 = diffmat(n, 2);
L = kron(D2, I) + kron(I, D2);
A = [L(ii,:) ; N + eta*II(bb,:)];
A = A(:,nc);
uu = u(xx,yy);
gg = (N + eta*II(bb,:)) * uu(:);
ff = f(xx,yy);
X = A \ [zeros(numIntPts, numBdyPts) ff(ii) ; eye(numBdyPts) zeros(numBdyPts, 1)];

% Interpolation operator for corner values:
Xii = xx(1,2:(n-1)).';
B = [Xii-1, -1-Xii].'; B(:,1:2:end) = -B(:,1:2:end);
if ( mod(n-1, 2) )
    B(2,:) = -B(2,:);
end
S = zeros(n^2, numBdyPts+1);
S(nc,:) = X;
S([1 n],:) = B*X(1:n-2,:);
S([end-n+1 end],:) = B*X(end-n+3:end,:);

G = (N - eta*II(bb,:));
Y = G * S(:,1:end-1);

norm(uu(:) - S*[gg(:) ; 1])
cond(Y)
figure(2), plot(eig(Y), 'o')
shg
