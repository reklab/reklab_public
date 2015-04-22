% State space Model Identification Toolbox.
% Version 1.0e    06-August-1999
%
% SMI functions
%
%   dordpo     - Compress data and estimate order for PO-MOESP.
%   dordpi     - Compress data and estimate order for PI-MOESP.
%   dordrs     - Compress data and estimate order for RS-MOESP.
%   dordeiv    - Compress data and estimate order for EIV-MOESP.
%   destac     - Estimate A and C.
%   destbdx    - Estimate B,D and X.
%   destk      - Estimate Kalman gain.
%   destx      - Estimate initial state.
%
% Optimization functions
%   dslslin    - Separable Least Squares for discrete time linear models.
%   dslswie    - Separable Least Squares for discrete time Wiener models.
%
% Utilities
%   vaf        - Variance Accounted For.
%   prbn       - Pseudo random binary noise.
%   shave      - Remove spikes from measured signals.
%   tchebest   - estimate MIMO static nonlinear function.
%   tchebsim   - simulate MIMO static nonlinear function.
%   ss2thon    - Translate state-space model to parameter vector.
%   th2sson    - Translate parameter vector to output normal state-space model.
%
% Demos
%   smidemo1   - Simple example with SISO system.
%   smidemo2   - Example with multiple batches and MIMO system.
%   smidemo3   - Example with SLS optimization.
%   smidemo4   - Example with Wiener model identification.
