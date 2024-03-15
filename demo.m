%% A Fast and Provable Algorithm for Sparse Phase Retrieval (MATLAB Package)

%% The algorithm consists of two stages:
% 1.Initialization: The first stage generates an initial estimate 
%                   using sparse spectral initialization method proposed in [1].
% 2.Refinement: The second stage refines the initial estimate to obtain the ground truth signal
%               using our proposed algorithm.

%% Reference 
% [1] Gauri Jagatap and Chinmay Hegde, "Sample-Efficient Algorithms for Recovering Structured Signals From Magnitude-Only Measurements," 
%     IEEE Transactions on Information Theory, vol. 65, no. 7, pp. 4434-4456, July 2019.

close all;
clc;
clear;
addpath('utils','functions')
rng(100)

%% problem setting
n = 10000;     % signal dimension
m = 5000;      % number of measurements
s = 100;       % sparsity of the signal

%% generate signal
z         = generate_true_signal(n,s);  % generate ground truth signal
[y_abs,A] = measure_signal(m,z);        % generate senging matrix and phaseless measurements

%% signal reconstruction
opts.maxiter  = 1e2;
opts.true     = z;
opts.toltrue  = 1e-6;
opts.display  = 1;
out           = solver_newton_projection(y_abs,A,s,opts);         

 

