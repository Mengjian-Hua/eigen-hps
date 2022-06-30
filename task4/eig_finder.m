function [V,D] = eig_finder(dom,pdo)
% This overloads Matlab built-in function eigs and it
% finds the eigenvalues (stored in a diagonal matrix D)
% and eigenvectors (columns of V) of the operator pdo defined
% in a closed domain dom


lambda = -1; % initial guess
previous_eig = lambda; % previously found eigenvalue

N_eigs = 10; % assume that the total number of eigenvalues is 3.
num_eigvalues_found = 0; % number of found eigenvalues
D = zeros(N_eigs,N_eigs); % matrix that stores all eigenvalues
V = zeros(N_eigs,N_eigs); % matrix that stores all eigenvectors
num_jumps = 0; % somtimes it may require more than one jump

while (num_eigvalues_found < N_eigs)
    
    lambda_new = root_finder(lambda,dom,pdo,1e-4);
    % if we get stuck in a piece, then try to go to the next piece
    if(abs(lambda_new - previous_eig)<1e-3)
        jump = 1*(num_jumps+1);
        lambda = lambda - jump; % try to go to another piece
        num_jumps = num_jumps + 1;
    else
        % if we enter into a new piece, then find the eigenvalue within it
        lambda = lambda_new; 
        num_eigvalues_found = num_eigvalues_found + 1;
        D(num_eigvalues_found,num_eigvalues_found) = lambda;
        previous_eig = lambda;
        num_jumps = 0; 
    end
   
end

end