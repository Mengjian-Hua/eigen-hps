function D = Dirichlet_diffmat(N,order)

    % we implement Dirichlet boundary conditions here 
    D = diffmat(N, order);
    D(1,:) = 0;
    D(1,1) = 1;
    D(end,:) = 0;
    D(end,end) = 1;
    
    
end