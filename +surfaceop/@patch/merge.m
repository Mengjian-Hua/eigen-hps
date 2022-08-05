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
[J1, J2, J31, J32, flip1, flip2, dom, edges] = intersect(a, b);

% Extract ItI maps:
Ra = a.R; Rb = b.R;
R33_a = Ra(J31,J31);
R33_b = Rb(J32,J32);
W = inv(1-R33_b*R33_a);
R_tau = [Ra(J1,J1) + Ra(J1,J31)*W*R33_b*Ra(J31,J1)...
        -Ra(J1,J31)*W*Rb(J32,J2);...
        -Rb(J2,J32)*(Ra(J31,J1)+R33_a*W*R33_b*Ra(J31,J1))...
        Rb(J2,J2) + Rb(J2,J32)*R33_a*W*Rb(J32,J2)];


z = [ scl2.*flip1*Ra(s1,i1) scl1.*flip2*Rb(s2,i2) ];
z_part = scl2.*flip1*a.du_part(s1,:) + scl1.*flip2*b.du_part(s2,:);

% Check for a closed surface at the top level
if ( rankdef && isempty(i1) && isempty(i2) )
    % Fix rank deficiency with Leslie's ones matrix trick:
    A = A + ones(size(A));
end

% Store the decomposition for reuse in updateRHS():
dA = decomposition(A);
S = dA \ z;
u_part = dA \ z_part;

if ( isIllConditioned(dA) )
    keyboard
end

% Compute new D2N maps:
Z12 = zeros(numel(i1), numel(i2));
%                                 |--- rhs ----|
% D2N = [ D2Na(i1,i1),  Z12,         D2Na(i1,end) ;
%         Z12.',        D2Nb(i2,i2), D2Nb(i2,end) ] ...
%     + [ D2Na(i1,s1)*flip1.' ; D2Nb(i2,s2)*flip2.' ] * S;
M = [ Ra(i1,s1)*flip1.' ; Rb(i2,s2)*flip2.' ];
D2N = [ Ra(i1,i1) Z12 ; Z12.' Rb(i2,i2) ] + M * S;
du_part = [ a.du_part(i1,:) ; b.du_part(i2,:) ] + M * u_part;

% Construct the new patch:
xyz = [a.xyz(i1,:) ; b.xyz(i2,:)];
id = [a.id ; b.id];
c = surfaceop.parent(dom, id, S, D2N, D2N_scl, u_part, du_part, A, dA, edges, xyz, a, b, ...
    {i1, s1}, {i2, s2}, flip1, flip2, scl1, scl2);

end
