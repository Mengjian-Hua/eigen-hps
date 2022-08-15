function u = solve(P, bc)
%SOLVE   Solve a leaf patch.
%   U = SOLVE(P, BC) returns a cell array containing the solution values U.

% Extract the domain from the patch:
id = P.id;
n = size(P.domain.x{id}, 1);

% if ( ~isnumeric(bc) )
%     % Evaluate the RHS if given a function handle:
%     bc = feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3));
% elseif ( isscalar(bc) )
%     % Convert a scalar to a constant vector:
%     bc = repmat(bc, size(P.xyz, 1), 1);
% end

% Evaluate the solution operator for the patch:
%u = P.S * [bc ; 1]; % The 1 accounts for the particular part.
u = P.X*[bc(:);1];
u = reshape(u, n, n);

% Return cell output for consistency with PARENT/SOLVE():
u = {u};

end
