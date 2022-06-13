function D = diffmat_periodic(N)

% This function creates a Chebyshev differentation matrix with periodic
% boundary conditions. (i.e. function values and first derivatives match
% at the two endpoints)
% Installation of Chebfun is required. 

D = diffmat(N,2);
D1 = D(1,:);
D(1,:) = 0;
D(1,1) = 1;
D(1,end) = -1;
D(end,:) = D1 - D(end,:);





end