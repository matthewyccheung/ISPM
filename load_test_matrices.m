%% Load Test Matrices

clc
close all
clear all

bcsstk15 = struct2sparsemat(load('bcsstk15.mat'));
mahindas = struct2sparsemat(load('mahindas.mat'));
nos3 = struct2sparsemat(load('nos3.mat'));
west0479 = struct2sparsemat(load('west0479.mat'));

figure;
subplot(2,2,1);
imshow(full(bcsstk15));
title("bcsstk15");
subplot(2,2,2);
imshow(full(mahindas));
title("mahindas");
subplot(2,2,3);
imshow(full(nos3));
title("nos3");
subplot(2,2,4);
imshow(full(west0479));
title("west0479");

function Y = struct2sparsemat(X)
    Y = X.Problem.A;
end
