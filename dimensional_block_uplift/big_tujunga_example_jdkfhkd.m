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
S = trunk(S);

% find the outlet of the stream 
Sout = streampoi(S,'outlets','ix');

% define the drainage basin and crop the DEM
DB = drainagebasins(FD,Sout);
DEM.Z(DB.Z == 0) = nan;
DEM = crop(DEM);

% invert streams with a K. Note that these numbers are arbitrary
K = 5e-6;
tau_inc = 2e5;
[A,U,S,Stau,tau_steps] = linear_inversion_block_uplift_erodibility(DEM,K,tau_inc);

Sz = DEM.Z(S.IXgrid);
Sz = sortrows(Sz)-min(Sz);
t = sortrows(Stau);
figure
subplot(2,1,1)
plot(t./1e6,A*U,'k-','linewidth',3); hold on
plot(t./1e6,Sz,'g-','linewidth',2)
legend('modeled','observed')
xlabel('Distance');
ylabel('Elevation')

Szs = movmean(Sz,100);
slo = (Szs(2:end)-Szs(1:end-1))./(t(2:end)-t(1:end-1));

subplot(2,1,2)
stairs(tau_steps./1e6,U.*1000,'color',[0 0 0],'lineWidth',3); hold on
plot(t(1:end-1)./1e6,slo.*1000,'g.')
xlabel('Distance'); ylabel('Slope');
legend('modeled','observed')