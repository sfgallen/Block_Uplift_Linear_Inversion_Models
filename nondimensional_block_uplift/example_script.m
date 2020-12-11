% example script for running inversion code
clear
clc

% add path to topotoolbox folder
addpath(genpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-2.2'));

%% Load a DEM as a topotoolbox GRIDobj
% load clipped watershed-scale DEM as a gridobj
DEM = load('DEM_ws_1.mat');
DEM = DEM.tDEM;

% invert streams
[A,Umod,S,Schi,chi_steps] = linear_inversion_block_uplift(DEM,'n_inc',10);

% example of how to convert from nondimentional to natural units and plot
K = 5e-6;
Ao = 1;
m = 0.45;

figure(99)
stairs(chi_steps./(K*Ao^m)./1e6,Umod.*(K*Ao^m).*1000);
xlabel('\tau (Myr)'); ylabel('Uplift rate (mm/yr)');