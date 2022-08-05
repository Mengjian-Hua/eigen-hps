classdef leaf < surfaceop.patch
%SURFACEOP.LEAF   Leaf subclass of a patch (where subproblems are solved).
%   P = SURFACEOP.LEAF(DOMAIN, S, D2N, EDGES, AINV) creates a LEAF object
%   P and assigns each of the inputs to their associated properties in P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTES:
% * Note: We assume that a LEAF has four sides and that its boundary nodes
%   XYZ are stored in the order "left", "right", "down", "up".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties

        n        % Discretization size of patch.
        % X        % Local (homogeneous BC) solution operator.
        normal_d % Normal derivative operator.
                 % (These are stored so the RHS can be efficiently updated)
    end

    methods

        function P = leaf(dom, id, X, R, u_part, Iu_part, edges, D2N, normal_d)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.n = size(dom.x{id}, 1);  % Discretization size.
            P.domain = dom;            % Domain.
            P.id = id;                 % Index of patch in domain.
            P.X = X;                   % Solution operator.
            P.R = R;                   % ItI map
            P.u_part = u_part;
            P.Iu_part = Iu_part;       % Outgoing impedance data (particular)
            P.edges = edges;           % Edges.
%             P.xyz = xyz;               % Boundary nodes.
            P.D2N = D2N;                  % Recovered D2N map
            P.normal_d = normal_d;     % Normal derivative operator.
            P.len = 1;

        end

    end

    methods ( Static )

        % Initialize an array of LEAF objects.
        P = initialize(op, dom, rhs);

    end

end
