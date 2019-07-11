%% Steepest Descent (SD) Method

clc
close all
clear all

run('load_test_matrices.m');
mats = ["bcsstk15", "mahindas", "nos3", "west0479"];
maxiter = 100; tol = eps;

figure; set(gca,'yscale', 'log'); hold on;
for i = 1:length(mats)
    A = eval(mats(i));
    n = size(A, 1);
    b = randi(100, n, 1);
    actual_sol = A\b;
    xk = ones(n, 1); 
    [xk, norm_rks] = sd(A, b, xk, maxiter, tol);
    iters = linspace(0, size(norm_rks, 2), size(norm_rks, 2));
    plot(iters, norm_rks);  
end
title('Steepest Descent (SD) Method Residuals Comparison for Test Matrices');
legend(mats);
xlabel('Iteration');
ylabel('Residual');
xlim([0, maxiter]);


function [xk, norm_rks] = sd(A, b, x0, maxiter, tol)
    iter = 0; xk = x0; x_prev = inf; norm_rks = [];
    while (norm(x_prev - xk) >= tol) & (iter <= maxiter)
        iter = iter + 1;
        rk = b - A*xk;
        norm_rks = [norm_rks, norm(rk)/(norm(A, 1)*norm(xk) + norm(b))];
        ak = rk'*rk/(rk'*A*rk);
        xprev = xk;
        xk = xk + ak*rk;
    end
end