% example script for running inversion code
clear
clc

% add path to topotoolbox folder
addpath(genpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-2.2'));

%% Load a DEM as a topotoolbox GRIDobj
% load DEM and clip to largest catchment for inversion as a gridobj
DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');

% fill sinks in DEM
DEM = fillsinks(DEM);

% make a flowobj
FD = FLOWobj(DEM);

% define stream network and largest connected component
S = STREAMobj(FD,'minarea',1e6/(DEM.cellsize)^2);
S = klargestconncomps(S,1);

% find the outlet of the stream 
Sout = streampoi(S,'outlets','ix');

% define the drainage basin and crop the DEM
DB = drainagebasins(FD,Sout);
DEM.Z(DB.Z == 0) = nan;
DEM = crop(DEM);


% invert streams with a K. Note that these numbers are arbitrary
K = 5e-6;
tau_inc = 1e5;
[A,U,S,Stau,tau_steps] = linear_inversion_block_uplift_erodibility(DEM,K,tau_inc);