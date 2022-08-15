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
        X        % Local (homogeneous BC) solution operator.
        normal_d % Normal derivative operator.
                 % (These are stored so the RHS can be efficiently updated)
        u_part
        du_part
    end

    methods

        function P = leaf(dom, id, R, D2N, D2N_scl, u_part, du_part, edges, xyz, w, X, normal_d,eta)

            % Construct empty patch:
            if ( nargin == 0 )
                return
            end

            % Assign domain and operators:
            P.n = size(dom.x{id}, 1); % Discretization size.
            P.domain = dom;           % Domain.
            P.id = id;                % Index of patch in domain.
            P.R = R;                  % ItI map
            P.D2N = D2N;              % Dirichlet-to-Neumann map.
            P.D2N_scl = D2N_scl;      % Scalings for Dirichlet-to-Neumann map.
            P.u_part = u_part;        % Particular solution.
            P.du_part = du_part;      % Impedance data of particular solution.
            P.edges = edges;          % Boundary edges.
            P.xyz = xyz;              % Boundary nodes.
            P.w = w;                  % Boundary quadrature weights.
            P.X = X;                  % Local solution operator.
            P.normal_d = normal_d;    % Normal derivative operator.
            P.len = 1;
            P.eta = eta;
        end

    end

    methods ( Static )

        % Initialize an array of LEAF objects.
        P = initialize(op, dom, rhs);

    end

end
