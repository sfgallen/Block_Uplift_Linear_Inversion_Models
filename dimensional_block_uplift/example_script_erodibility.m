% example script for running inversion code
clear
clc

% add path to topotoolbox folder
addpath(genpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-master'));

%% Load a DEM as a topotoolbox GRIDobj
% load clipped watershed-scale DEM as a gridobj
DEM = load('DEM_ws_1.mat');
DEM = DEM.tDEM;

% invert streams with a K
K = 5e-6;
tau_inc = 1e5;
[A,U,S,Stau,tau_steps] = linear_inversion_block_uplift_erodibility(DEM,K,tau_inc);
