%% Task 1 (using Chebfun)
clear all;close all;clc
N = 1e3; % number of Chebyshev nodes
K = 20; % number of eigenvalues we want to compare
% Homogeneous Dirichilet boundary condition
index = 1;
for N = 22:200
    D = diffmat(N, 2,'dirichlet', 'dirichlet');
    eigvalues = sort(eig(D)',"descend");
    k = 1:K;
    true_eig = -k.^2*pi^2/4;
    L1_error(index) = sum(abs(eigvalues(3:K+2)-true_eig))/(N-2);
    L2_error(index) = sqrt(sum((eigvalues(3:K+2)-true_eig).^2)/(N-2));
    L_inf_error(index) = max(abs(eigvalues(3:K+2)-true_eig));
    index = index+1;
end
% Periodic boundary condition
index = 1;
for N = 22:200
    D = diffmat(N, 2,'periodic');
    eigvalues = sort(eig(D)',"descend");
    k = 1:K;
    true_eig = -k.^2*pi^2;
    L1_error_periodic(index) = sum(abs(eigvalues(2:K+1)-true_eig))/(N-2);
    L2_error_periodic(index) = sqrt(sum((eigvalues(2:K+1)-true_eig).^2)/(N-2));
    L_inf_error_periodic(index) = max(abs(eigvalues(2:K+1)-true_eig));
    index = index+1;
end

%% Plotting the results
figure(1)
semilogy(22:200,L1_error,"linewidth",2);
hold on
semilogy(22:200,L2_error,"linewidth",2);
semilogy(22:200,L_inf_error,"linewidth",2);
set(gca,"fontsize",20)
legend({"$L^1$","$L^2$","$L^\infty$"},"interpreter","latex")
xlabel("Number of Chebyshev nodes")
ylabel("Error in different norms")
title("Homogeneous Dirichlet BCs (first 20 eigenvalues)")

