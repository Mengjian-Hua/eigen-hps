N = 20;
K = 10;
[A, B] = diffmat_glue(N, K, 2);

% Generalized eigenvalue problem:
e = eigs(A, B, 10, 'sm');
sqrt(abs(e / pi^2 * 4))

%% Reorder degrees of freedom to get a block diagonal matrix in the top left
bb = 1:(N-1):K*N-1;            % Glue indices
ii = 1:K*(N-1)+1; ii(bb) = []; % Interior indices
AA = [ A(ii,ii) A(ii,bb) ;
       A(bb,ii) A(bb,bb) ];
BB = [ B(ii,ii) B(ii,bb) ;
       B(bb,ii) B(bb,bb) ];
spy(AA)

% As we only shuffled rows and columns, this has the same eigenvalues:
e = eigs(AA, BB, 10, 'sm');
sqrt(abs(e / pi^2 * 4))

%% Take Schur complement to eliminate continuity conditions
S = A(ii,ii) - A(ii,bb) * (A(bb,bb) \ A(bb,ii));

% Now B = I, so we have a standard eigenvalue problem:
e = eigs(S, 10, 'sm');
sqrt(abs(e / pi^2 * 4))
