function x = null_overload(A)

%% Solve for the singular linear system Ax = b
%  based on the randomized method in the paper by
%  Sifuentes, Gimbutas, & Greengard
%  https://arxiv.org/pdf/1401.3068.pdf

dim = size(A);
p = randn(dim(1),1);
q = randn(dim(1),1);

x = randn(dim(1),1);
b = A*x;
y = (A + p*q')\b;
x = x-y;


end