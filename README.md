# bldp - Bregman Log Determinant Preconditioner

This package implements the preconditioners in [1] in MATLAB.

[1] [Bock, Andreas A., Martin S. Andersen. "A New Matrix Truncation Method for Improving Approximate Factorisation Preconditioners." arXiv preprint arXiv:2312.06417 (2023).](https://arxiv.org/abs/2312.06417)

## Installation

The examples we show depend on SuiteSparse. Download and place it in the
directory root.

## Repository structure

.
├── bldp.m              # Main library for computing the preconditioners from [1].
│                       # Implements exact truncations and approximations using 
│                       # Nyström and Krylov-Schur.
├── example_bldp.m      # Contains two small examples of `bldp.m`
├── run_pcg_small.m     # Run experiments from [1] for smaller matrices (exact truncations)
├── run_pcg_small.m     # Run experiments from [1] for smaller matrices (uses Krylov-Schur and Nyström to approximate the truncations)
├── run_pcg_large.m     # Run experiments from [1] for large matrices (uses Krylov-Schur and Nyström to approximate the truncations)
├── utils               # SuiteSparse helper functions and various plotting functions
├── tests               # Contains some tests
├── LICENSE
└── README.md           # This file