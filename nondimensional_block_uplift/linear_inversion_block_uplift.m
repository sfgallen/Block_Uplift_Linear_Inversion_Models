function [A,Umod,S,Schi] = linear_inversion_block_uplift(DEM,varargin)
% Block uplift inversion code after Goren et al. (2014) JGR-Earth Surface.
% Inputs:
%       (1) DEM clipped to a watershed as a TopoToolbox GRIDobj
%       (2) crita - critical drainage area for channel head initiation in
%           m^2 (default = 1e6)
%       (3) mn - the m over n ratio (concavity) of river system (default =
%           0.5)
%       (4) Ao - reference drainage area for chi-integration (default =
%           1)
%       (5) n_inc - number of time increments solved for in inversion
%           (default = 10)
%       (6) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes that DEM is already 
%          filled or carved}
%
% Outputs:
%       (1) A - forward model matrix
%       (2) Umod - recovered uplift history
%       (3) S - topotoolbox STREAMobj
%       (4) Schi - chi of the stream network
%
% Author: Sean F. Gallen
% Date Modified: 03/10/2020
% email: sean.gallen[at]colostate.edu

%%Parse Inputs
p = inputParser;         
p.FunctionName = 'linear_inversion_block_uplift';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));

% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'Ao', 1, @(x) isscalar(x));
addOptional(p,'n_inc', 10, @(x) isscalar(x));
addOptional(p,'flowOption', []);

parse(p,DEM, varargin{:});
DEM   = p.Results.DEM;
crita = p.Results.crita;
mn = p.Results.mn;
Ao = p.Results.Ao;
n_inc = p.Results.n_inc;

%% Process grids
% set nan values if it hasn't already been done
DEM.Z(DEM.Z <= -9999) = NaN;

% declare cellsize
cs = DEM.cellsize;

% make direction FLOWobj
% flow routing options
if isempty(p.Results.flowOption)
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'fill')
    DEM = fillsinks(DEM);
    FD = FLOWobj(DEM);
elseif strcmp(p.Results.flowOption, 'carve')
    FD = FLOWobj(DEM,'preprocess','carve');
    DEM = imposemin(FD,DEM);
else
    error('fillOption is not "fill" or "carve"');
end

% make flow accumulation grid
DA = flowacc(FD).*(FD.cellsize^2);

% Create stream network object (S)
S  = STREAMobj(FD,'minarea',crita/(cs)^2);
S = klargestconncomps(S,1);

Sz = DEM.Z(S.IXgrid);

Schi = MakeChiNetwork(S,DA,mn,Ao);

%% coding Goren et al., 2014 block uplift inversion model
% sort by elevation
table = sortrows([Sz,Schi]);
z_vec = table(:,1)-min(table(:,1));
chi_vec = table(:,2);

n = length(z_vec);  % number of nodes

% define desired time increments.
min_chi = floor(min(chi_vec));
max_chi = ceil(max(chi_vec));

% make time vector
inc_chi = (max_chi-min_chi)/n_inc;
chi_steps = min_chi:inc_chi:max_chi; % time steps
q = length(chi_steps);

% load the time matrix with data
A = zeros(n,q);
for i = 2:q
    inc_chi_t=chi_steps(i)-chi_steps(i-1);
    A(chi_vec >= chi_steps(i),i-1) = inc_chi_t;
end
for i = 1:n
    loc_sum = sum(A(i,:));
    loc = find(A(i,:)==0,1);
    A(i,loc) = chi_vec(i) - loc_sum; %this is the reminder
end

% calculate Upri
ZoverAsum = z_vec./sum(A,2);
ZoverAsum(ZoverAsum == Inf) = nan;
Upri = ones(q,1)*(1/n)*nansum(ZoverAsum);

% make q by q identity matrix
I = eye(q,q);

% impose dampening (e.g. smoothing) parameter, (Gamma)
Gam = logspace(-2,5,1000);

%% find the best-fit Gamma and plot results
cols = jet(length(Gam));
MisFit = nan(1,length(Gam));

figure()
subplot(2,2,1)
for i = 1:length(Gam)
    % now we can invert to find U following Goren from Tarantola, 1987
    Umod = Upri + (A'*A + Gam(i)^2*I)\A'*(z_vec-A*Upri);
%     subplot(2,1,1);
%     stairs(chi_steps,Umod,'color',cols(i,:),'lineWidth',0.5); hold on
    MisFit(i) = (1/(n-q))*sqrt(sum((z_vec - A*Umod).^2));
    %subplot(2,1,2)
    plot(1./Gam(i),MisFit(i),'.','color',cols(i,:)); hold on
end

[ind,~] = turingPointFinder(1./Gam,MisFit);
%subplot(2,1,2)
plot(1/Gam(ind),MisFit(ind),'ko');
xlabel('1/\Gamma'); ylabel('normalized misfit');

% define best-fit Gamma and do the final inversion
Gam = Gam(ind);

%% invert for uplift history with best damping factor
Umod = Upri + (A'*A + Gam^2*I)\A'*(z_vec-A*Upri);

% plot uplift results
subplot(2,2,2)
% incremental uplift
stairs(chi_steps,Umod,'color',[0 0 0],'lineWidth',2); hold on
xlabel('\chi'); ylabel('U*');

% cummulative uplift
upWRTtime = cumtrapz(chi_steps,Umod);
subplot(2,2,4)
stairs(chi_steps,upWRTtime,'color',[0 0 0],'lineWidth',1); hold on
xlabel('\chi'); ylabel('total baselevel fall(m)');

% plot observed and best-fit model results
subplot(2,2,3)
plot(chi_vec,z_vec,'.','color',[0.5 0.5 0.5]); hold on
plot(chi_vec,A*Umod,'k.');%,'color',Pcols(k,:));
xlabel('\chi'); ylabel('elevation (m)')
legend({'observed','modeled'});
end

function [idx,perpDist] = turingPointFinder(x,y)
% function finds "corner" on a plot with x and y asymptotes
%
% Author: Sean F. Gallen
% Date modified: 03/23/2016
% email: sean.gallen{at}colostate.edu

% Two endpoints on the curve "data"
xend = [x(1) x(end)];
yend = [y(1) y(end)];
% The slope of the line connecting the two endpoints
m = ( yend(2) - yend(1) )/( xend(2) - xend(1) );
% Point on the curve (xc,yc), point on the line (xl,yl)
perpDist = zeros(length(x),1);
for i = 1:length(x)
    xc = x(i) ; yc = y(i);
    yl = ( (m * xc) + (m^2 * yc) - (m * xend(1)) + yend(1) )/(1+ m^2);
    xl = xc - m*(yl - yc);
    % distance^2
    d2 = (xl - xc)^2 + (yl - yc)^2;
    perpDist(i) = (d2);
end
[~, idx] = max(perpDist);
end

function [Schi] = MakeChiNetwork(S,DA,mn,Ao)
%
% MakeChiNetwork.m will make calculate chi for a TopoToolbox stream
% network.
%
% Inputs:
% S         - TopoToolbox STREAMobj.
% DA        - Drainage area grid IN MAP UNITs (e.g. m^2) as a GRIDobj.
% mn        - m/n or refence concavity (theta) value.
%
% Outputs:
% Schi      - Chi for the stream network.
%
% Author: Sean F. Gallen
% Date modified: 11/19/2018
% email: sean.gallen{at}colostate.edu

p = inputParser;
p.FunctionName = 'MakeChiNetwork';
addRequired(p,'S', @(x) isa(x,'STREAMobj'));
addRequired(p,'A', @(x) isa(x,'GRIDobj'));
addRequired(p,'mn', @(x) isscalar(x));
addRequired(p,'Ao', @(x) isscalar(x));

parse(p,S,DA,mn,Ao);

% get variables ready for chi integration
Schi = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = (Ao./(DA.Z(S.IXgrid))).^mn;       % chi transformation variable

h = waitbar(0,'calculating \chi for full stream network...');
% calculating chi and tau_star for the entire river network
for lp = numel(Six):-1:1
    Schi(Six(lp)) = Schi(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end
close(h);
end