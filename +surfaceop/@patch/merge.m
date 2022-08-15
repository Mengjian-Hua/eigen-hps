function c = merge(a, b, rankdef)
%MERGELR   Merge two patch objects.
%   C = MERGE(A, B) returns a patch C formed by merging the two patches A
%   and B. Typically A and B will be adjacent and any common edges will be
%   eliminated by enforcing continuity and continuity of the derivative
%   across the boundary.

% Parse inputs:
if ( nargin == 0 )
    c = [];
    return
elseif ( nargin == 1 )
    c = a;
    return
elseif ( nargin == 2 )
    rankdef = false;
end

% Compute the indices of intersecting points in a and b.
[J1, J2, J31, J32,flip1, flip2, scl1, scl2, D2N_scl, dom, edges] = intersect(a, b);

% Extract ItI maps:
Ra = a.R; Rb = b.R;
eta = a.eta;
R33_a = Ra(J31,J31);
R33_b = Rb(J32,J32);
W = inv(eye(size(R33_b*R33_a))-R33_b*R33_a);
R = [Ra(J1,J1) + Ra(J1,J31)*W*R33_b*Ra(J31,J1)...
        -Ra(J1,J31)*W*Rb(J32,J2);...
        -Rb(J2,J32)*(Ra(J31,J1)+R33_a*W*R33_b*Ra(J31,J1))...
        Rb(J2,J2) + Rb(J2,J32)*R33_a*W*Rb(J32,J2)];
D2N = -eta*(R-eye(size(R)))\(R+eye(size(R)));
S_alpha = [W*Rb(J32,J32)*Ra(J31,J1) -W*Rb(J32,J2)];
S_beta = -[Ra(J31,J1)+Ra(J31,J31)*W*Rb(J32,J32)*Ra(J31,J1) W*Rb(J32,J2)];

% Extract D2N maps:
% D2Na = a.D2N; D2Nb = b.D2N;
% Compute new solution operator:
% - The Dirichlet-to-Neumann maps on singular elements have their Jacobians
%   factored out. Therefore, we need to multiply the continuity conditions
%   by the multiplication matrices SCL1 and SCL2. The coordinate maps are
%   such that multiplying D2NA/B by SCL1/2 cancels out, so D2NA/B should
%   only be multiplied by SCL2/1, respectively.
% A = -( scl2.*flip1*D2Na(J31,J31)*flip1.' + scl1.*flip2*D2Nb(J32,J32)*flip2.' );
% z = [ scl2.*flip1*D2Na(J31,J1) scl1.*flip2*D2Nb(J32,J2) ];
% z_part = scl2.*flip1*a.du_part(J31,:) + scl1.*flip2*b.du_part(J32,:);
% 
% % Check for a closed surface at the top level:
% if ( rankdef && isempty(J1) && isempty(J2) )
%     % Fix rank deficiency with Leslie's ones matrix trick:
%     w = a.w(J31);
%     A = A + w*w'; % or is it sqrt(w)*sqrt(w)' ?
% end
% 
% % Store the decomposition for reuse in updateRHS():
% dA = decomposition(A);
% S = dA \ z;
% u_part = dA \ z_part;
% % 
% % % Compute new D2N maps:
% Z12 = zeros(numel(J1), numel(J2));
% M = [ D2Na(J1,J31)*flip1.' ; D2Nb(J2,J32)*flip2.' ];
% D2N = [ D2Na(J1,J1) Z12 ; Z12.' D2Nb(J2,J2) ] + M * S;
% du_part = [ a.du_part(J1,:) ; b.du_part(J2,:) ] + M * u_part;

% Construct the new patch:
xyz = [a.xyz(J1,:) ; b.xyz(J2,:)];
w = [a.w(J1) ; b.w(J2)];
id = [a.id ; b.id];
c = surfaceop.parent(dom, id, R, D2N, D2N_scl, S_alpha,S_beta, ...
    edges, xyz, w, a, b, {J1, J31}, {J2, J32}, flip1, flip2, scl1, scl2);

end
