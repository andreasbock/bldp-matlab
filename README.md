# bldp - Bregman Log Determinant Preconditioner

This package implements the preconditioners in [1] in MATLAB.

[1] [Bock, Andreas A., Martin S. Andersen. "A New Matrix Truncation Method for Improving Approximate Factorisation Preconditioners." arXiv preprint arXiv:2312.06417 (2023).](https://arxiv.org/abs/2312.06417)

## Installation

The examples we show depend on SuiteSparse. Download and place it in the
directory root.

## Repository structure

.
├── example_truncations.m      # Contains a small example illustrating the main idea in presented in [1].
├── example_bldp.m             # Solve a few linear systems using the various preconditioners from `bldp.m`. 
├── bldp.m                     # Main library for computing the preconditioners from [1]. Implements
│                              # exact truncations and approximations using randomised Nyström and Krylov-Schur.
├── run_pcg_small.m            # Run experiments from [1] for smaller matrices (compares exact truncations
│                              # to Krylov-Schur/Nyström approximations).
├── run_pcg_large.m            # Run experiments from [1] for large matrices (uses Krylov-Schur and Nyström
│                              # to approximate the truncations).
├── utils                      # SuiteSparse helper functions and various plotting functions.
├── tests                      # Contains some tests.
├── LICENSE
└── README.md                  # This file.