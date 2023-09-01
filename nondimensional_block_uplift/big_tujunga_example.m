% example script for running inversion code
clear
clc

% add path to topotoolbox folder
addpath(genpath('C:\Users\sfgallen\Documents\topo_toolbox\topotoolbox-master'));

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

% invert streams
[A,Umod,S,Schi,chi_steps,res] = linear_inversion_block_uplift(DEM,'n_inc',10,'mn',0.45);

% example of how to convert from nondimentional to natural units and plot
% Note that the K value is arbitrary
K = 5e-6;
Ao = 1;
mn = 0.5;

figure(99)
stairs(chi_steps./(K*Ao^mn)./1e6,Umod.*(K*Ao^mn).*1000);
xlabel('\tau (Myr)'); ylabel('Uplift rate (mm/yr)');


