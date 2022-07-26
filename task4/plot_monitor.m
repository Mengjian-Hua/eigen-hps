% Plot the monitor function for the Laplace-Beltrami operator

clear all;
l = 5;
n = 20;
dom = surfacemesh.sphere(n, 2);
pdo = [];
pdo.lap = 1;
n = 100;
lambda = linspace(-1,-13,n);


for i = 1:n
    
    f(i) = eval_monitor(lambda(i),pdo,dom);
    
end

figure(1)
plot(lambda,f,"linewidth",2)
set(gca,"fontsize",15)
title("Monitor function for the Laplace-Beltrami operator on a sphere")
xlabel("$\lambda$","interpreter","latex")
ylabel("Monitor function")
xlim([-13 -1])
hold on;plot(lambda,zeros(size(lambda)),"linewidth",2)