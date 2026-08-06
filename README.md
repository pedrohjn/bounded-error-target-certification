# Finite-Time Target Certification under Bounded-Error Position Sensing

This repository contains the MATLAB code and numerical results
associated with the paper:

> Pedro H. J. Nardelli, Pedro E. Gória Silva, and Nicola Marchetti,
> "Finite-Time Target Certification in Multi-Agent Systems under
> Bounded-Error Position Sensing."

The code evaluates finite-time certification in multi-agent
multi-target systems under bounded-error position sensing.

## Main results reproduced

The repository reproduces the numerical results concerning:

1. pairwise certification probability;
2. system-level certification probability;
3. conditional normalized certification time;
4. early-certification probability;
5. heterogeneous per-agent sensing-error bounds;
6. different error distributions within a common deterministic bound.

## Requirements

- MATLAB R2021b or later is recommended.
- No additional MATLAB toolboxes are required.

The code relies on standard MATLAB functions for random-number
generation, tables, plotting, and file export.

## Repository structure

```text
code/
    run_numerical_analysis.m
    plot_numerical_results.m
    run_sensor_extensions.m
    sample_disk_error.m

results/
    numerical_results.csv
    sensor_extensions_results.csv

figures/
    numerical_tradeoffs.pdf
    sensor_extensions.pdf
