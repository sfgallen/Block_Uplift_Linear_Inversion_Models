% example script for running inversion code
clear
clc

% add path to topotoolbox folder
addpath(genpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-2.2'));

%% Load a DEM as a topotoolbox GRIDobj
% load clipped watershed-scale DEM as a gridobj
DEM = load('DEM_ws_1.mat');
DEM = DEM.tDEM;

% invert streams with a K
K = 1e-5;
tau_inc = 6e4;
[A,U,S,Stau] = linear_inversion_block_uplift_erodibility(DEM,K,tau_inc);
