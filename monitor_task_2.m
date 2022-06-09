%% Task 2

clear all;close all;clc
N = 1e3; % number of Chebyshev nodes

% Homogeneous Dirichilet boundary condition
A = diffmat(N, 2,'dirichlet', 'dirichlet');
A_lambda = @(lambda) A - lambda*eye(N);
w = randn(N,1); % random vector
v = randn(N,1); % random vector
F1 = @(lambda) 1./w'*(A_lambda(lambda)\v); % First monitor function
K = 10; % number of random vectors
wk = randn(N,K); % a set of random vectors
vk = randn(N,K); % a set of random vectors
F2 = @(lambda) 1./trace((wk'*(A_lambda(lambda)\vk)).^2);% second monitor function
F3 = @(lambda) 1./sum((A_lambda(lambda)\vk).^2/N,"all");% third monitor function
%% evaluting the monitor functions
lambda = -4:0.01:-1;
for i = 1:length(lambda)
    
    F1_lambda(i) = F1(lambda(i)); 
    F2_lambda(i) = F2(lambda(i)); 
    F3_lambda(i) = F3(lambda(i)); 
    
end
%% Plotting the results
figure(1)
plot(lambda,F1_lambda,"linewidth",2)
title("The first monitor function")
figure(2)
plot(lambda,F2_lambda,"linewidth",2)
title("The second monitor function")
figure(3)
plot(lambda,F3_lambda,"linewidth",2)
title("The third monitor function")

%% root finding with the bisection method for the second monitor function
%  which is the smoothest one among all three monitor functions

% initial guess (i.e. the eigenvalue lies within interval [-4,-1])
guess_1 = -4;
guess_2 = -1;

tol = 1e-3; % tolerance

left = F2(guess_1);
right = F2(guess_2);

num_evaluation = 2; % number of evalution of the monitor function

while(abs(guess_1-guess_2)>tol)
    
    if(abs(left) > abs(right))
        guess_1 = (guess_1 + guess_2)/2;
        left = F2(guess_1);
    else
        guess_2 = (guess_1 + guess_2)/2;
        left = F2(guess_2);
    end
    num_evaluation = num_evaluation + 1;
end

eigvalue = (guess_1 + guess_2)/2;

%% root finding with the fixed-point iteration for the second monitor function

F = @(lambda) F2(lambda) - lambda;
guess = -1;
diff = inf; % lambda_{n+1} - lambda_n
tol = 1e-3;

num_evaluation = 0;

while (diff>tol)
    
    guess = F(guess);
    num_evaluation = num_evaluation + 1;
    
end

eigvalue = guess;
