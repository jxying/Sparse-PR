# SparsePR

This MATLAB package provides implementations of a fast second-order algorithm based on Newton projection [1] for **sparse phase retrieval**.

The algorithm consists of two stages:

1. Initialization: The first stage generates an initial estimate that is close to the ground truth signal, using the sparse spectral initialization method proposed in [2].

2. Refinement: The second stage refines the initial estimate to obtain the ground truth signal, using our proposed second-order algorithm.


## To run the code:
(1) Download the source files.
(2) Run 'demo.m' in MATLAB.


## References
[1] Jian-Feng Cai, Yu Long, Ruixue Wen, and Jiaxi Ying, "A Fast and Provable Algorithm for Sparse Phase Retrieval", in nternational Conference on Learning Representations (ICLR), 2024.

[2] Gauri Jagatap, and Chinmay Hegde, "Sample-Efficient Algorithms for Recovering Structured Signals From Magnitude-Only Measurements," in IEEE Transactions on Information Theory, vol. 65, no. 7, pp. 4434-4456, July 2019.
