
function P = cheb_projection(m,n,kind)
% This function returns the m*n matrix 
% which interpolates data from a grid 
% of n*4 Chebyshev nodes of the first kind to a5r
% grid of 4*m-4 Chebyshev nodes of the second kind
% with every one of the four corners being averaged from 
% the two edges that share it.

if ( nargin == 0 )
    test();
    return
end

if(kind == "second") % convert the first kind nodes to the second kind nodes
    [x, ~, w] = chebpts(n, 1);
    y = chebpts(m, 2);
    B = barymat(y, x, w); % Barycentric interpolation matrix for one edge

    P = zeros(4*m-4, 4*n);
    P(1:m,1:n)                   = B;
    P(m:2*m-1,n+1:2*n)           = B;
    P(2*m-1:3*m-2,2*n+1:3*n)     = B;
    P([3*m-2:4*m-4 1],3*n+1:end) = B;

    % Average the corners
    P(1,:) = P(1,:)/2;
    P(m,:) = P(m,:)/2;
    P(2*m-1,:) = P(2*m-1,:)/2;
    P(3*m-2,:) = P(3*m-2,:)/2;
elseif(kind == "first")
    [x, ~, w] = chebpts(n, 2);
    y = chebpts(m, 1);
    B = barymat(y, x, w); % Barycentric interpolation matrix for one edge
    
    % P1,P2,P3,P4 are interpolation matrices for each edge. 
    P1 = zeros(4*m, 4*n-4); 
    P2 = zeros(4*m, 4*n-4);
    P3 = zeros(4*m, 4*n-4);
    P4 = zeros(4*m, 4*n-4);
    P1(1:m,1:n) = B;
    P2(m+1:2*m,n:2*n-1) = B;
    P3(2*m+1:3*m,2*n-1:3*n-2) = B;
    P4(3*m+1:4*m,[3*n-2:end 1]) = B;
    P = P1 + P2 + P3 + P4;
elseif(kind == "first_duplicate")
    [x, ~, w] = chebpts(n, 2);
    y = chebpts(m, 1);
    B = barymat(y, x, w); % Barycentric interpolation matrix for one edge
    P = blkdiag(B,B,B,B);

end

end

function test()

n = 30;

% N-1 first-kind Chebyshev nodes on the boundary, ordered counterclockwise
% from the top left
x1 = chebpts(n-1, 1);
y1 = chebpts(n-1, 1);
xx1 = [-1+0*x1 ; x1 ; 1+0*x1 ; flipud(x1)];
yy1 = [flipud(y1) ; -1+0*y1 ; y1 ; 1+0*y1];

% N second-kind Chebyshev nodes on the boundary, ordered counterclockwise
% from the top left
x2 = chebpts(n, 2);
y2 = chebpts(n, 2);
xx2 = [-1+0*x2 ; x2 ; 1+0*x2 ; flipud(x2)];
yy2 = [flipud(y2) ; -1+0*y2 ; y2 ; 1+0*y2];
xx2([n 2*n 3*n 4*n]) = [];
yy2([n 2*n 3*n 4*n]) = [];

% Projection matrix from first-kind skeleton to second-kind skeleton
P = cheb_projection(n, n-1,2);

u = randnfun2;
norm(u(xx2,yy2) - P*u(xx1,yy1))

% test the other direction
P = cheb_projection(n-1, n,1);
u = randnfun2;
norm(u(xx1,yy1) - P*u(xx2,yy2))

end
