function L = initialize(op, dom, rhs)
%INITIALIZE   Initialize an array of LEAF objects.
%   L = SURFACEOP.LEAF.INITIALIZE(OP, DOM) returns a cell array L of LEAF
%   objects which contain the solution and D2N operators for Poisson's
%   equation on the domain DOM with zero righthand side.
%
%   L = SURFACEOP.LEAF.INITIALIZE(OP, DOM, RHS) is as above, but with the
%   righthand side RHS, which may be a scalar or a function handle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( isempty(dom) )
    L = [];
    return
end

assert(isa(dom, 'surfacemesh'), 'Invalid domain.');

if ( nargin < 3 )
    % Default to homogeneous problem:
    rhs = 0;
end

numPatches = length(dom);
n = size(dom.x{1}, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% DEFINE REFERENCE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%

[X, Y] = chebpts2(n);             % Chebyshev points and grid.
ii = abs(X) < 1 & abs(Y) < 1;     % Interior indices.
ee = ~ii;
% ii([1,end],[1,end]) = true;     % Treat corners as interior.
eta = 2i;
% Note that here the index sets are different from what we had before
kindFrom = 1;
kindTo   = 1;
ss = [1:n 3*n-3:4*n-4 n:2:3*n-3 4*n-4 1 n+1:2:3*n-3];


[x0, ~, v0] = chebpts(n,   kindFrom);
[x1, ~, v1] = chebpts(n,   2);
[x2, ~, v2] = chebpts(n-2, kindTo);

nbdy = n;
numBdyPts = 4*nbdy; 
numIntPts = sum(ii(:));

xn1 = chebpts(nbdy, 1);
xxb = [ -1+0*xn1 ; 1+0*xn1 ; xn1 ; xn1 ];
yyb = [ xn1 ; xn1 ; -1+0*xn1 ; 1+0*xn1 ];

B = zeros(numBdyPts, n^2);
for k = 1:numBdyPts
    B(k,:) = kron(barymat(xxb(k), x0, v0), barymat(yyb(k), x0, v0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = diffmat(n);

I = eye(n);
II = kron(I, I);
Du = kron(D, I);
Dv = kron(I, D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANT RHS? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define scalar RHSs:
if ( isnumeric(rhs) && isscalar(rhs) )
    rhs_eval_full = repmat(rhs,n^2, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% SOLVE LOCAL PROBLEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Each patch needs a different solution operator.

% Initialize
L = cell(numPatches, 1);

% Loop over each patch:
for k = 1:numPatches

    % Define the left and right edges for this patch:
    x = dom.x{k};
    y = dom.y{k};
    z = dom.z{k};
    edges = [ x(1,1) y(1,1) z(1,1) x(n,1) y(n,1) z(n,1) n ;  % "Left" side
              x(1,n) y(1,n) z(1,n) x(n,n) y(n,n) z(n,n) n ;  % "Right" side
              x(1,1) y(1,1) z(1,1) x(1,n) y(1,n) z(1,n) n ;  % "Down" side
              x(n,1) y(n,1) z(n,1) x(n,n) y(n,n) z(n,n) n ]; % "Up" side

    % Evaluate non-constant RHSs if required:
    if ( isa(rhs, 'function_handle') )
        rhs_eval_full = feval(rhs, x, y, z);
    elseif ( isa(rhs, 'surfacefun') )
        rhs_eval_full = rhs.vals{k};
    elseif ( iscell(rhs) )
        rhs_eval_full = rhs{k};
    end

    %[Dx, Dy, Dz] = diffs(x, y, z);
    Dx = dom.ux{k}(:).*Du + dom.vx{k}(:).*Dv;
    Dy = dom.uy{k}(:).*Du + dom.vy{k}(:).*Dv;
    Dz = dom.uz{k}(:).*Du + dom.vz{k}(:).*Dv;
        
    Dxx = Dx^2;
    Dyy = Dy^2;
    Dzz = Dz^2;
    Dxy = Dx*Dy;
    Dyz = Dy*Dz;
    Dxz = Dx*Dz;
    if ( dom.singular(k) )
        J = dom.J{k}(:);
        Jx = Dx2*J; Jy = Dy*J; Jz = Dz*J;
        Dxx = J.*Dxx - Jx.*Dx;
        Dyy = J.*Dyy - Jy.*Dy;
        Dzz = J.*Dzz - Jz.*Dz;
        Dxy = J.*Dxy - Jx.*Dy;
        Dyz = J.*Dyz - Jy.*Dz;
        Dxz = J.*Dxz - Jz.*Dx;
        Dx = J.^2.*Dx;
        Dy = J.^2.*Dy;
        Dz = J.^2.*Dz;
        II = J.^3.*II;
    end

    for name = fieldnames(op).'
        name = name{1};
        if ( isa(op.(name), 'function_handle') )
            op.(name) = feval(op.(name), x, y, z);
        elseif ( isa(op.(name), 'surfacefun') )
            op.(name) = op.(name).vals{k};
        end
    end

    A = zeros(n^2,n^2);
    if ( op.dxx ~= 0 ), A = A + op.dxx(:).*Dxx; end
    if ( op.dyy ~= 0 ), A = A + op.dyy(:).*Dyy; end
    if ( op.dzz ~= 0 ), A = A + op.dzz(:).*Dzz; end
    if ( op.dxy ~= 0 ), A = A + op.dxy(:).*Dxy; end
    if ( op.dyz ~= 0 ), A = A + op.dyz(:).*Dyz; end
    if ( op.dxz ~= 0 ), A = A + op.dxz(:).*Dxz; end
    if ( op.dx  ~= 0 ), A = A + op.dx(:).*Dx;   end
    if ( op.dy  ~= 0 ), A = A + op.dy(:).*Dy;   end
    if ( op.dz  ~= 0 ), A = A + op.dz(:).*Dz;   end
    if ( op.b   ~= 0 ), A = A + op.b(:).*II;    end
    
    % construct the solution operator (ItI)
    [nl, nr, nd, nu] = normals(x, y, z);
    for i = 1:3
          
        nl(:,i) = chebvals2chebvals(nl(:,i), 2, 1);
        nr(:,i) = chebvals2chebvals(nr(:,i), 2, 1);
        nd(:,i) = chebvals2chebvals(nd(:,i), 2, 1);
        nu(:,i) = chebvals2chebvals(nu(:,i), 2, 1);

    end
    
    normal_d = zeros(numBdyPts, n^2);
    S02 = barymat(x2, x0, v0);

    S01 = barymat(x0, x1, v1);
    P01 = kron(S01, S01);
    

    S_rhs = barymat(x2, x1, v1);
    P_rhs = kron(S_rhs, S_rhs);
    rhs_eval = P_rhs * rhs_eval_full(:);
    
    normal_d(1:nbdy,:)  = nl(:,1).*(B(1:nbdy,:)*P01*Dx)  + nl(:,2).*(B(1:nbdy,:)*P01*Dy)  + nl(:,3).*(B(1:nbdy,:)*P01*Dz);
    normal_d(nbdy+1:2*nbdy,:) = nr(:,1).*(B(nbdy+1:2*nbdy,:)*P01*Dx) + nr(:,2).*(B(nbdy+1:2*nbdy,:)*P01*Dy) + nr(:,3).*(B(nbdy+1:2*nbdy,:)*P01*Dz);
    normal_d(2*nbdy+1:3*nbdy,:)  = nd(:,1).*(B(2*nbdy+1:3*nbdy,:)*P01*Dx)  + nd(:,2).*(B(2*nbdy+1:3*nbdy,:)*P01*Dy)  + nd(:,3).*(B(2*nbdy+1:3*nbdy,:)*P01*Dz);
    normal_d(3*nbdy+1:4*nbdy,:)    = nu(:,1).*(B(3*nbdy+1:4*nbdy,:)*P01*Dx)    + nu(:,2).*(B(3*nbdy+1:4*nbdy,:)*P01*Dy)  + nu(:,3).*(B(3*nbdy+1:4*nbdy,:)*P01*Dz);
    
    % Construct solution operator with impedance data:

    if ( dom.singular(k) )
        % solution operator, denoted by X, with incoming impedance data
        F = normal_d + eta*B*P01; % equation (2.9) 
        A = P_rhs*A;
        B1 = [A ; F];
        dB1 = decomposition(B1, 'cod');
        rhsX = [zeros(numIntPts, numBdyPts) rhs_eval(:); eye(numBdyPts) zeros(numBdyPts, 1)];
        X = dB1\rhsX; % equation below (2.10)
    else
        % solution operator, denoted by X, with incoming impedance data
        F = normal_d + eta*B*P01; % equation (2.9) 
        %A = P01*A;
        %B1 = [A(ii,:);F];
        A = P_rhs*A;
        B1 = [A ; F];
        dB1 = decomposition(B1);
        rhsX = [zeros(numIntPts, numBdyPts) rhs_eval(:); eye(numBdyPts) zeros(numBdyPts, 1)];
        X = dB1\rhsX; % equation below (2.10)
    end
      
    % Construct the ItI map
    G = normal_d - eta*B*P01;
    R = G*X(:,1:numBdyPts);
    [xn1, ~, vn1] = chebpts(nbdy, 1);
    [xn2, ~, vn2] = chebpts(n-2, 1);
    C  = barymat(xn2, xn1, vn1);
    C1 = barymat(xn1, xn2, vn2);
    CC = blkdiag(C, C, C, C);
    CC1 = blkdiag(C1, C1, C1, C1);
    R = CC * R * CC1;
    
    u_part = X(:,end);
    Iu_part = CC*G*u_part;

    if ( dom.singular(k) )
        D2N_scl = cell(4, 1);
        J = dom.J{k}.^3; Jee = J(ee); Jss = Jee(ss);
        D2N_scl{1} = Jss(leftIdx); % Left
        D2N_scl{2} = Jss(rightIdx);  % Right
        D2N_scl{3} = Jss(downIdx); % Down
        D2N_scl{4} = Jss(upIdx);  % Up
    else
        D2N_scl = {ones(n-2,1), ones(n-2,1), ones(n-2,1), ones(n-2,1)}.';
    end
    
    % Assemble the patch:
    xee = x(ee);
    yee = y(ee);
    zee = z(ee);
    xyz = [xee(ss) yee(ss) zee(ss)];

    for side = 1:4
        xyz(side*n-n+1:side*n,1) = chebvals2chebvals(xyz(side*n-n+1:side*n,1),2,1);
        xyz(side*n-n+1:side*n,2) = chebvals2chebvals(xyz(side*n-n+1:side*n,2),2,1);
        xyz(side*n-n+1:side*n,3) = chebvals2chebvals(xyz(side*n-n+1:side*n,3),2,1);
    end

    w = chebtech2.quadwts(n); w = w(:);
    ww = w .* w.' .* sqrt(dom.J{k});
    ww = ww(ee);
    ww = ww(ss);

    D2N = -eta*(R-eye(4*n-8))\(R+eye(4*n-8));
    L{k} = surfaceop.leaf(dom, k, R, D2N, D2N_scl, u_part, Iu_part, edges, xyz, ww, X, normal_d,eta);
    
    % test
%     u_true = spherefun.sphharm(1,0);
%     uu2 = u_true(x,y,z); % second kind nodes
%     gg = F*uu2(:); % compute the incoming impedance data
%     norm(uu2(:) - X*[gg(:);1])
%     gg = CC*gg; 
%     ff = CC*G*uu2(:); % outgoing impedance data 
%     norm(R*gg + CC*G*u_part - ff)
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% DEFINE OPERATORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Dx, Dy, Dz] = diffs(x, y, z)
%DIFFS   Differential operators on surfaces.
%   [DX, DY, DZ] = DIFFS(X, Y, Z) returns the tangential derivative
%   operators DX, DY, and DZ over the surface defined by the coordinates
%   (X, Y, Z). The coordinates must be given as N x N matrices sampled at
%   tensor-product Chebyshev nodes. DX, DY, and DZ are operators of size
%   N^2 x N^2.

persistent D
n = size(x, 1);
if ( size(D, 1) ~= n )
    D = diffmat(n);
end

xu = x * D.'; xv = D * x;
yu = y * D.'; yv = D * y;
zu = z * D.'; zv = D * z;
E = xu.*xu + yu.*yu + zu.*zu;
G = xv.*xv + yv.*yv + zv.*zv;
F = xu.*xv + yu.*yv + zu.*zv;
J = E.*G - F.^2;
ux = (G.*xu-F.*xv)./J; vx = (E.*xv-F.*xu)./J;
uy = (G.*yu-F.*yv)./J; vy = (E.*yv-F.*yu)./J;
uz = (G.*zu-F.*zv)./J; vz = (E.*zv-F.*zu)./J;
I = eye(n); % or speye?
Du = kron(D, I);
Dv = kron(I, D);
Dx = ux(:).*Du + vx(:).*Dv;
Dy = uy(:).*Du + vy(:).*Dv;
Dz = uz(:).*Du + vz(:).*Dv;

end

function [nl, nr, nd, nu] = normals(x, y, z)
%NORMALS   Outward pointing normal vectors to the edges of a mapping.

persistent D
n = size(x, 1);
if ( size(D, 1) ~= n )
    D = diffmat(n);
end

xu = x * D.'; xv = D * x;
yu = y * D.'; yv = D * y;
zu = z * D.'; zv = D * z;

nl = -normalize([xu(:,1)   yu(:,1)   zu(:,1)]);
nr =  normalize([xu(:,n)   yu(:,n)   zu(:,n)]);
nd = -normalize([xv(1,:).' yv(1,:).' zv(1,:).']);
nu =  normalize([xv(n,:).' yv(n,:).' zv(n,:).']);

tangent = normalize([xv(:,1) yv(:,1) zv(:,1)]);
nl = nl - tangent .* dot(nl, tangent, 2);
tangent = normalize([xv(:,n) yv(:,n) zv(:,n)]);
nr = nr - tangent .* dot(nr, tangent, 2);
tangent = normalize([xu(1,:).' yu(1,:).' zu(1,:).']);
nd = nd - tangent .* dot(nd, tangent, 2);
tangent = normalize([xu(n,:).' yu(n,:).' zu(n,:).']);
nu = nu - tangent .* dot(nu, tangent, 2);

nl = normalize(nl);
nr = normalize(nr);
nd = normalize(nd);
nu = normalize(nu);

% nl([1,n],:) = [];
% nr([1,n],:) = [];
% nd([1,n],:) = [];
% nu([1,n],:) = [];

end

function v = normalize(v)

v = v ./ sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);

end
