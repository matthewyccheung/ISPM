%% Conjugate Gradient (CG) Method

clc
close all
clear all

run('load_test_matrices.m');
A = bcsstk15;
n = size(A, 1);
b = randi(10, n, 1);
x0 = b/norm(b);



run('load_test_matrices.m');
mats = ["bcsstk15", "mahindas", "nos3", "west0479"];
maxiter = 100; tol = eps;

figure; set(gca,'yscale', 'log'); hold on;
for i = 1:length(mats)
    A = eval(mats(i));
    n = size(A, 1);
    b = randi(100, n, 1);
    actual_sol = A\b;
    x0 = b/norm(b);
    [x, norm_rks] = conjgrad(A, b, x0, maxiter);
    iters = linspace(0, size(norm_rks, 2), size(norm_rks, 2));
    plot(iters, norm_rks);  
end
title('Conjugate Gradient (CG) Method Residuals Comparison for Test Matrices');
legend(mats);
xlabel('Iteration');
ylabel('Residual');
xlim([0, 100]);
hold off;



function [x, norm_rks] = conjgrad(A, b, x, maxiter)
    r = b - A*x;
    d = r;
    r2_prev = r'*r;
    norm_rks = [r2_prev];
    for i = 1:maxiter
        Ad = A*d;
        alpha = r2_prev/(d'*Ad);
        r = r - alpha*Ad;
        x = x + alpha*d;
        norm_rks = [norm_rks, norm(r)/(norm(A, 1)*norm(x) + norm(b))];
        r2 = r' * r;
        if r2 < 1e-10
              break;
        end
        beta = r2/r2_prev;
        d = r + beta*d;
        r2_prev = r2;
    end
end