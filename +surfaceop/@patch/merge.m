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
W_inv = eye(size(R33_b*R33_a))-R33_b*R33_a;
R = [Ra(J1,J1) + Ra(J1,J31)*W*R33_b*Ra(J31,J1)...
        -Ra(J1,J31)*W*Rb(J32,J2);...
        -Rb(J2,J32)*(Ra(J31,J1)+R33_a*W*R33_b*Ra(J31,J1))...
        Rb(J2,J2) + Rb(J2,J32)*R33_a*W*Rb(J32,J2)];
D2N = -eta*(R-eye(size(R)))\(R+eye(size(R)));
S_alpha = [W*Rb(J32,J32)*Ra(J31,J1) -W*Rb(J32,J2)];
S_beta = -[Ra(J31,J1)+Ra(J31,J31)*W*Rb(J32,J32)*Ra(J31,J1) W*Rb(J32,J2)];
Iu_part_a = a.Iu_part;
Iu_part_b = b.Iu_part;

if(isempty(J31) || isempty(J1))
    Iu_part = [Iu_part_a;Iu_part_b];
else
    Iu_part = [Iu_part_a(J1);Iu_part_b(J2)] + [Ra(J1,J31)*W_inv*(R33_a*Iu_part_a(J31) - Iu_part_b(J32));...
            -Rb(J2,J32)*(eye(size(R33_a*W_inv*R33_b)) + R33_a*W_inv*R33_b)*Iu_part_a(J31)...
            + Rb(J2,J32)*R33_a*W_inv*Iu_part_b(J32)];
end

% Construct the new patch:
xyz = [a.xyz(J1,:) ; b.xyz(J2,:)];
w = [a.w(J1) ; b.w(J2)];
id = [a.id ; b.id];
% if(length(J1) == length(J31))
%     keyboard
% end
c = surfaceop.parent(dom, id, R, D2N, D2N_scl, Iu_part,S_alpha,S_beta, ...
    edges, xyz, w, a, b, {J1, J31}, {J2, J32}, flip1, flip2, scl1, scl2);

end
