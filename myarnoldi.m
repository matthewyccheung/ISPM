%% Arnoldi Iteration

clc
close all
clear all

run('load_test_matrices.m');
A = west0479;
n = size(A, 1);
v = ones(n, 1);
m = 30;
[Q1, H1] = arnoldi_orthog(A, v, m);
[Q2, H2] = arnoldi(A, v, m);


% 3
norm(A*Q1(:, 1:end-1) - Q1*H1)
norm(A*Q2(:, 1:end-1) - Q2*H2)
norm(eye(size(Q1, 2), size(Q1,2)) - Q1'*Q1)
norm(eye(size(Q2, 2), size(Q2,2)) - Q2'*Q2)

% 4
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
figure; hold on;
scatter(real(eigenvalues_H1), imag(eigenvalues_H1), 'r+');
scatter(real(eigenvalues_A), imag(eigenvalues_A));
title('Hessenberg vs Exact Eigenvalues of A for iterations = 30');
xlabel('Real'); ylabel('Imag');
legend('Ritz', 'Exact');
norm(A*Q1(:, 1:end-1) - Q1*H1)

% 5
m = 10;
[Q1, H1] = arnoldi_orthog(A, v, m);
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:10);
figure; hold on;
scatter(real(eigenvalues_H1), imag(eigenvalues_H1), 'r+');
scatter(real(eigenvalues_A), imag(eigenvalues_A));
title('Hessenberg vs Exact Eigenvalues of A for iterations = 10');
xlabel('Real'); ylabel('Imag');
legend('Ritz', 'Exact');

m = 20;
[Q1, H1] = arnoldi_orthog(A, v, m);
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:20);
figure; hold on;
scatter(real(eigenvalues_H1), imag(eigenvalues_H1), 'r+');
scatter(real(eigenvalues_A), imag(eigenvalues_A));
title('Hessenberg vs Exact Eigenvalues of A for iterations = 20');
xlabel('Real'); ylabel('Imag');
legend('Ritz', 'Exact');

%6
m = 30;
[Q1, H1] = arnoldi(A, v, m);
figure; hold on;
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for iterations = 30')
legend('H1', 'A'); hold off;
norm(eigenvalues_A - eigenvalues_H1)

m = 60;
[Q1, H1] = arnoldi(A, v, m);
figure; hold on;
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:60);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for iterations = 60')
legend('H1', 'A'); hold off;
norm(eigenvalues_A - eigenvalues_H1)



m = 30;
[Q1, H1] = arnoldi_orthog(A, v, m);
figure; hold on;
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for iterations = 30')
legend('H1', 'A'); hold off;
% norm(eye(size(Q1, 2), size(Q1,2)) - Q1'*Q1)
norm(eigenvalues_A - eigenvalues_H1)

m = 60;
[Q1, H1] = arnoldi_orthog(A, v, m);
figure; hold on;
eigenvalues_H1 = eig(H1(1:end-1, :));
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:60);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for iterations = 60')
legend('H1', 'A'); hold off;
% norm(eye(size(Q1, 2), size(Q1,2)) - Q1'*Q1)
norm(eigenvalues_A - eigenvalues_H1)



% 7
tau = 10;
m = 30;
[Q, H, U, mu] = arnoldi_shift_invert(A, v, m, tau);
figure; hold on;
eigenvalues_H1 = mu;
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for tao = 10 and iterations = 30');
legend('H1', 'A'); hold off;

tau = -10;
m = 30;
[Q, H, U, mu] = arnoldi_shift_invert(A, v, m, tau);
figure; hold on;
eigenvalues_H1 = mu;
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for tao = -10 and iterations = 30');
legend('H1', 'A'); hold off;

tau = -5;
m = 30;
[Q, H] = arnoldi_shift_invert(A, v, m, tau);
figure; hold on;
eigenvalues_H1 = 1./eig(H(1:end-1, :)) + tau;
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for tao = -5 and iterations = 30');
legend('H1', 'A'); hold off;

tau = 5;
m = 30;
[Q, H, U, mu] = arnoldi_shift_invert(A, v, m, tau);
figure; hold on;
eigenvalues_H1 = 1./eig(H(1:end-1, :)) + tau;
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for tao = 5 and iterations = 30');
legend('H1', 'A'); hold off;

tau = 1;
m = 30;
[Q, H] = arnoldi_shift_invert(A, v, m, tau);
figure; hold on;
eigenvalues_H1 = 1./eig(H(1:end-1, :)) + tau;
eigenvalues_A = eig(full(A)); 
eigenvalues_A = eigenvalues_A(1:30);
scatter(real(eigenvalues_A), imag(eigenvalues_A), 'r+');
scatter(real(eigenvalues_H1), imag(eigenvalues_H1));
title('Eigenvalues of A for tao = 1 and iterations = 30');
legend('H1', 'A'); hold off;



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


function [Q, H] = arnoldi_orthog(A, v, m)
    n = size(A, 1);
    H = zeros(m+1, m);
    Q = zeros(n, m+1);
    Q(:, 1) = v/norm(v);
    for k = 1:m
        v = A*Q(:, k);
        
        for j = 1:k
            H(j, k) = Q(:, j)'*v;
            v = v - H(j, k)*Q(:, j);
            v_norm_prev = norm(v);
        end
 
        for j=1:k
           tmp = Q(:,j)'*v;
           v = v - tmp * Q(:, j);
           H(j, k) = H(j, k) + tmp;
        end
        
        H(k+1, k) = norm(v);
        Q(:, k+1) = v/H(k+1, k);
        
        if H(k+1, k) <= ((n*eps)*v_norm_prev)
            return;
        end
    end
end


function [Q, H, U, mu] = arnoldi_shift_invert(A, v, m, tau)
    n = size(A, 1);
    H = zeros(m+1, m);
    Q = zeros(n, m+1);
    Q(:, 1) = v/norm(v);
    U = zeros(n, 1);
    mu = zeros(n, 1);
    
    [L, R, P] = lu(A - tau*speye(n));
    
    for k = 1:m
        v = inv(P'*L*R)*Q(:, k);
        
        U(:, k) = v/norm(v);
        mu(k) = U(:, k)'*A*U(:, k);
        
        for j = 1:k
            H(j, k) = Q(:, j)'*v;
            v = v - H(j, k)*Q(:, j);
            v_norm_prev = norm(v);
        end
 
        for j=1:k
           tmp = Q(:,j)'*v;
           v = v - tmp * Q(:, j);
           H(j, k) = H(j, k) + tmp;
        end
        
        H(k+1, k) = norm(v);
        Q(:, k+1) = v/H(k+1, k);
        
        
        if H(k+1, k) <= ((n*eps)*v_norm_prev)
            return;
        end
    end
end
