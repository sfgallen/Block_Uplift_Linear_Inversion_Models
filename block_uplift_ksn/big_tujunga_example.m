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
Ksn_Inv = linear_inversion_ksn(DEM,'chi_inc',1,'mn',0.45,'gam',10);



% plot uplift results
subplot(3,1,1)
% incremental uplift
stairs(Ksn_Inv.chi_steps,Ksn_Inv.ksn_mod,'color',[0 0 0],'lineWidth',2); hold on
xlabel('\chi'); ylabel('ksn');

% plot observed and best-fit model results
subplot(3,1,2)
plot(Ksn_Inv.Schi,Ksn_Inv.Sz,'.','color',[0.5 0.5 0.5]); hold on
plot(Ksn_Inv.Schi,Ksn_Inv.A*Ksn_Inv.ksn_mod,'k.');%,'color',Pcols(k,:));
xlabel('\chi'); ylabel('elevation (m)')
legend({'observed','modeled'});

subplot(3,1,3)
plot(Ksn_Inv.S.distance./1000,Ksn_Inv.Sz,'.','color',[0.5 0.5 0.5]); hold on
plot(Ksn_Inv.S.distance./1000,Ksn_Inv.A*Ksn_Inv.ksn_mod,'k.');%,'color',Pcols(k,:));
xlabel('Distance (km)'); ylabel('elevation (m)')
legend({'observed','modeled'});


