function [D,B] = diffmat_glue(N,K,order)

% This function creates a Chebyshev differentation matrix with periodic
% boundary conditions. (i.e. function values and first derivatives match
% at the two endpoints)
% Installation of Chebfun is required. 

D = zeros(N*K-K+1,N*K-K+1);
B = eye(size(D));

    for i = 1:K
        
        % differentiation matrices
        first_index = 1+(i-1)*(N-1);
        last_index = first_index + (N-1);
        xi = -1+(first_index-1)*2/K;
        D(first_index:last_index,first_index:last_index) = diffmat(N,order,[xi xi+2/K]);
        
        % matching the first derivatives
        if(i>1)
            D1_1 = diffmat(N,1,[xi xi+2/K]);
            D1_2 = diffmat(N,1,[xi+2/K xi+4/K]);
            D1_first = D1_1(1,:);
            D1_last = D1_2(end,:);
            D(first_index,first_index+1:last_index) = -D1_first(1:end-1);
            D(first_index,first_index - (N-1):first_index-1) = D1_last(2:end);
            D(first_index,first_index) = D1_last(1) - D1_first(end);
        end
        
        
        B(first_index,:) = 0;
        B(last_index,:) = 0;
    end
    
%% boundary conditions at two endpoints
D(1,:) = 0;
D(1,1) = 1;

D(end,:) = 0;
D(end,end) = 1;
    


end
