classdef parent < surfaceop.patch
%SURFACEOP.PARENT   Parent subclass of a patch.
%   P = SURFACEOP.PARENT(DOMAIN, S, D2N, EDGES, XY, CHILD1, CHILD2, IDX1,
%   IDX2, L2G1, L2G2) creates a SURFACEOP.PARENT object P and assigns each
%   of the inputs to their associated properties in P.

    properties

        child1 = [] % Child patch
        child2 = [] % Child patch
        idx1        % How p.xyz relates to p.child1.xyz
        idx2        % How p.xyz relates to p.child2.xyz
        flip1
        flip2
        scl1
        scl2
        S1           % Impedance data operator for the first child
        S2           % Impedance data operator for the second child

    end

    methods

        function P = parent(domain, id, R, D2N, D2N_scl, Iu_part, S1, S2, edges, xyz, w, child1, child2, idx1, idx2, flip1, flip2, scl1, scl2)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.domain = domain;
            P.id = id;
            P.R = R;
            P.D2N = D2N;
            P.D2N_scl = D2N_scl;
%             P.u_part = u_part;
            P.Iu_part = Iu_part;
            P.S1 = S1;
            P.S2 = S2;
            P.edges = edges;
            P.xyz = xyz;
            P.w = w;
            P.len = child1.len + child2.len;
            P.eta = child1.eta;

            % Assign children:
            P.child1 = child1;
            P.child2 = child2;
            P.idx1 = idx1;
            P.idx2 = idx2;
            P.flip1 = flip1;
            P.flip2 = flip2;
            P.scl1 = scl1;
            P.scl2 = scl2;
            
        end

    end

end