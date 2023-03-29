function u = solve(P, bc)
%SOLVE   Solve a leaf patch.
%   U = SOLVE(P, BC) returns a cell array containing the solution values U.

% Extract the domain from the patch:
id = P.id;
n = size(P.domain.x{id}, 1);

eta = 2i;
if ( ~isnumeric(bc) && ~isa(bc,"surfacefun"))
    % Evaluate the RHS if given a function handle:
    bc = eta*feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3)) - abs(eta)^2*P.D2N*feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3));
elseif ( isscalar(bc) )
    % Convert a scalar to a constant vector:
    bc = eta*repmat(bc, size(P.xyz, 1), 1) + P.D2N*repmat(bc, size(P.xyz, 1), 1);
end

% Evaluate the solution operator for the patch:
%u = P.S * [bc ; 1]; % The 1 accounts for the particular part.
% [xn1, ~, vn1] = chebpts(n, 2);
% [xn2, ~, vn2] = chebpts(n, 1);
% C  = barymat(xn1, xn2, vn2);
% CC = blkdiag(C, C, C, C);


u = P.X*[bc(:);1];
u = reshape(u, n, n);

% Return cell output for consistency with PARENT/SOLVE():
u = {u};

end
