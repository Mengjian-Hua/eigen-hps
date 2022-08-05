classdef ( Abstract ) patch
%SURFACEOP.PATCH   Abstract patch object.

    properties

        domain  % Domain of patch.
        id      % Index of patch in domain.
        X       % Solution operator for patch.
        R     % Dirichlet-to-Neumann map for patch.
        u_part
        Iu_part
        edges   % Boundary edges of patch.
        % xyz     % Boundary grid points of patch.
        len
        D2N

    end

    methods ( Abstract )

        % Solve a patch.
        [u, d] = solve(P, bc);

        % Update RHS of a patch.
        f = updateRHS(f, rhs);

        % Number of degrees of freedom in a patch.
        N = numel(P);

        % Number of patches in a patch.
        N = length(P);

    end

end
