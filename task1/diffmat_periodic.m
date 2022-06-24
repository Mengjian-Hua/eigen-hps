function [D,B] = diffmat_periodic(N,order)

% This function creates a Chebyshev differentation matrix with periodic
% boundary conditions. (i.e. function values and first derivatives match
% at the two endpoints)
% Installation of Chebfun is required. 

D = diffmat(N,order);
D1 = diffmat(N,1);
D(1,:) = 0;
D(1,1) = 1;
D(1,end) = -1;
D(end,:) = D1(1,:) - D1(end,:);

B = zeros(N,N);
B(2:end-1,2:end-1) = eye(N-2);


end

%% test code
% [D,B] = diffmat_periodic(100,2);
% lambda = flip(sort(eig(D,B)));