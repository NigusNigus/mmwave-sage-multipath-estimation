 ## Beamsteering-Based SAGE for mmWave Multipath Estimation

This repository contains MATLAB code for estimating multipath channel
parameters in mmWave systems using a beamsteering-based SAGE algorithm.

The implementation works on beamformed (beamspace) measurements obtained
from coarse beamscanning and estimates AoD, AoA, delay, and complex path
gains.

The implementation follows the same formulation and equations as the
referenced paper, with minor implementation-level updates for clarity
and consistency.

## Reference
Wireless Multipath Component Estimation for Millimeter-Wave Beam-Scanning Systems 
IEEE, 2024.  
https://ieeexplore.ieee.org/abstract/document/10683217

 ## Files
- `SAGE_mmWave_Multipath_Estimation.m` – main script
- `SAGE.m` – SAGE loop
- `Construction_ofMPC.m` – beamsteering and delay-domain model
- `ExpectationXL.m` – E-step computation
- `make_initial_guess.m` – Parameter initialization
- `Z_Corr_thetaTX.m` – AoD update
- `Z_Corr_thetaRX.m` – AoA update
- `Z_Corr_Delay.m` – delay update
- `plot_final_results.m` – plotting results

 ## To Run

SAGE_mmWave_Multipath_Estimation
