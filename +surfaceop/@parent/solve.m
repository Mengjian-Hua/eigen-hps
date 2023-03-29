function u = solve(P, bc)
%SOLVE   Solve a parent patch.
%   U = SOLVE(P, BC) returns a cell array U of function values representing
%   the PDE solution on the parent P with Dirichlet boundary data given by
%   BC.
eta = 2i;
if ( ~isnumeric(bc) && ~isa(bc,"surfacefun"))
    % Evaluate the RHS if given a function handle:
    bc = eta*feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3)) - abs(eta)^2*P.D2N*feval(bc, P.xyz(:,1), P.xyz(:,2), P.xyz(:,3));

elseif ( isscalar(bc) )
    % Convert a scalar to a constant vector:
    bc = eta*repmat(bc, size(P.xyz, 1), 1) + P.D2N*repmat(bc, size(P.xyz, 1), 1);
end

% Evaluate the solution operator for the parent:
% u = P.S * [bc ; 1]; % The 1 accounts for the particular part.
Iu = conj(P.Iu_part); % convert outgoing ItI data to incoming ItI data


% Construct boundary conditions for children and recurse.

% Construct boundary indices:
i1 = cat(1,P.idx1{1});
i31 = cat(1,P.idx1{2});
i2 = cat(1,P.idx2{1});
i32 = cat(1,P.idx2{2});

bci1 = 1:numel(P.idx1{1});
bci2 = 1:numel(P.idx2{1});
if ( ~isempty(i1) )
    bci2 = bci2 + bci1(end);
end
% CAT() is 10x faster than CELL2MAT().
idx1 = cat(1, P.idx1{:}); % idx1 = cell2mat(P.idx1.');
idx2 = cat(1, P.idx2{:}); % idx2 = cell2mat(P.idx2.');

% Assemble boundary conditions for child patches:
if (isempty(P.R)) % top level
    ubc1 = zeros(size(P.child1.S2, 2), 1);
    ubc2 = zeros(size(P.child2.S2, 2), 1);
    if(isempty(bc))
        ubc1(idx1) = Iu(1:end/2);
        ubc2(idx2) = Iu(end/2+1:end);
    else
        ubc1(i1) = bc(bci1) - Iu(1:end/2);
        ubc2(i2) = bc(bci2) - Iu(end/2+1:end);
        ubc1(i31) = P.flip1.'*P.S1*(bc-Iu);
        ubc2(i32) = P.flip2.'*P.S2*(bc-Iu);
    end
elseif (isempty(P.S1) && isempty(P.S2))
        ubc1 = zeros(size(P.child1.S2, 2), 1);
        ubc2 = zeros(size(P.child1.S2, 2), 1);
        ubc1(idx1) = bc(bci1) - Iu(1:end/2);
        ubc2(idx2) = bc(bci2) - Iu(end/2+1:end);
else
        ubc1 = zeros(size(P.child1.xyz,1), 1);
        ubc2 = zeros(size(P.child1.xyz,1), 1);
        ubc1(i1) = bc(bci1) - Iu(1:end/2);
        ubc2(i2) = bc(bci2) - Iu(end/2+1:end);
        ubc1(i31) = P.flip1.'*P.S1*(bc-Iu);
        ubc2(i32) = P.flip2.'*P.S2*(bc-Iu);
end


% Solve for the child patches:
u1 = solve(P.child1, ubc1);
u2 = solve(P.child2, ubc2);

% Concatenate for output:
u = [u1 ; u2]; 

end
