%% Restarted GMRES Method

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
    [xk, norm_rks] = gm(A, b, xk, maxiter);
    iters = linspace(0, size(norm_rks, 2), size(norm_rks, 2));
    plot(iters, norm_rks);  
end
title('Restarted GMRES Method Comparison for Test Matrices');
legend(mats);
xlabel('Iteration');
ylabel('Residual');
xlim([0, maxiter]);

function [xk, norm_rks] = gm(A, b, xk, maxiter)
    norm_rks = [];
    n = maxiter;
    for i = 1:100
        [Q, H] = arnoldi(A, xk, maxiter);
        be1 = zeros(n+1, 1); be1(1) = norm(b);
        y = H\be1;
        xk = xk + Q(:, 1:end-1)*y;
        norm_rks = [norm_rks, norm(H*y - be1)/norm(b)];
    end
end


function [Q, H] = arnoldi(A, v, m)
    n = size(A, 1);
    H = zeros(m+1, m);
    Q = zeros(n, m+1);
    Q(:, 1) = v/norm(v);
    for k = 1:m
        v = A*Q(:, k);
        for j = 1:k
            H(j, k) = Q(:, j)'*v;
            v = v - H(j, k)*Q(:, j);
        end
        H(k+1, k) = norm(v);
        Q(:, k+1) = v/H(k+1, k);
    end
end