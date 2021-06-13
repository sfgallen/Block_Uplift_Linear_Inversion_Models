function [A,Umod,S,Stau,tau_steps] = linear_inversion_block_uplift_erodibility(DEM,K, tau_inc,varargin)
% Block uplift inversion code after Goren et al. (2014) JGR-Earth Surface.
% Inputs:
%       Required:
%       (1) DEM clipped to a watershed as a TopoToolbox GRIDobj
%       (2) Erodibility constant, K, calculated assuming n = 1
%       (3) Tau increment in years (divides network into even tau incs)
%
%       Optional:
%       (3) crita - critical drainage area for channel head initiation in
%           m^2 (default = 1e6)
%       (3) mn - the m over n ratio (concavity) of river system (default =
%           0.5)
%       (6) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes that DEM is already 
%          filled or carved}
%
% Outputs:
%       (1) A - forward model matrix
%       (2) Umod - recovered uplift history in m/yr
%       (3) S - topotoolbox STREAMobj
%       (4) Stau - tau of the stream network in yrs
%       (5) tau_steps - increments of chi used in inversion in yrs
%
% Author: Sean F. Gallen
% Date Modified: 12/11/2020
% email: sean.gallen[at]colostate.edu

%%Parse Inputs
p = inputParser;         
p.FunctionName = 'linear_inversion_block_uplift_erodibility';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));
addRequired(p,'K', @(x) isscalar(x));
addRequired(p,'tau_inc', @(x) isscalar(x));


% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'flowOption', []);

parse(p,DEM, K, tau_inc, varargin{:});
DEM   = p.Results.DEM;
K = p.Results.K;
tau_inc = p.Results.tau_inc;
crita = p.Results.crita;
mn = p.Results.mn;

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
%S = trunk(S);

Sz = DEM.Z(S.IXgrid);

Stau = MakeTauNetwork(S,K,DA,mn);

%% coding Goren et al., 2014 block uplift inversion model
% sort by elevation
table = sortrows([Sz,Stau]);
z_vec = table(:,1)-min(table(:,1));
tau_vec = table(:,2);

n = length(z_vec);  % number of nodes

% define desired time increments.
min_tau = floor(min(tau_vec));
max_tau = ceil(max(tau_vec));

if tau_inc > max_tau
    error(['Error: Your Tau increment is larger than the maximum Tau. ',...
        'The maximum Tau is: ', num2str(max_tau),' years. Reduce the the increment',...
        'so that you have ~5 to 10 tau increments between 0 and the maximum Tau.']);
end

% make time vector
tau_steps = min_tau:tau_inc:max_tau; % time steps
q = length(tau_steps);

% load the time matrix with data
A = zeros(n,q);
for i = 2:q
    inc_chi_t=tau_steps(i)-tau_steps(i-1);
    A(tau_vec >= tau_steps(i),i-1) = inc_chi_t;
end
for i = 1:n
    loc_sum = sum(A(i,:));
    loc = find(A(i,:)==0,1);
    A(i,loc) = tau_vec(i) - loc_sum; %this is the reminder
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
    
    MisFit(i) = (1/(n-q))*sqrt(sum((z_vec - A*Umod).^2));
    %subplot(2,1,2)
    plot(1./Gam(i),MisFit(i),'.','color',cols(i,:)); hold on
end

[ind,~] = turingPointFinder(1./Gam,MisFit);

plot(1/Gam(ind),MisFit(ind),'ko','markerfacecolor',[0.5 0.5 0.5]);
xlabel('1/\Gamma'); ylabel('normalized misfit');

% define best-fit Gamma and do the final inversion
Gam = Gam(ind);

%% invert for uplift history with best damping factor
Umod = Upri + (A'*A + Gam^2*I)\A'*(z_vec-A*Upri);

% plot uplift results
subplot(2,2,2)
% incremental uplift
stairs(tau_steps./1e6,Umod.*1000,'color',[0 0 0],'lineWidth',2); hold on
xlabel('\tau (Myr)'); ylabel('U (mm/yr)');

% cummulative uplift
upWRTtime = cumtrapz(tau_steps,Umod);
subplot(2,2,4)
stairs(tau_steps./1e6,upWRTtime,'color',[0 0 0],'lineWidth',1); hold on
xlabel('\tau (Myr)'); ylabel('total baselevel fall(m)');

% plot observed and best-fit model results
subplot(2,2,3)
plot(tau_vec./1e6,z_vec,'.','color',[0.5 0.5 0.5]); hold on
plot(tau_vec./1e6,A*Umod,'k.');%,'color',Pcols(k,:));
xlabel('\tau (Myr)'); ylabel('elevation (m)')
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

function [Stau] = MakeTauNetwork(S,K,DA,mn)
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
addRequired(p,'K', @(x) isscalar(x));
addRequired(p,'DA', @(x) isa(x,'GRIDobj'));
addRequired(p,'mn', @(x) isscalar(x));


parse(p,S,K,DA,mn);

% get variables ready for chi integration
Stau = zeros(size(S.distance));
Six = S.ix;
Sixc = S.ixc;
Sx = S.distance;
Sa = (1./(DA.Z(S.IXgrid))).^mn;       % chi transformation variable

h = waitbar(0,'calculating \tau for full stream network...');
% calculating chi and tau_star for the entire river network
for lp = numel(Six):-1:1
    Stau(Six(lp)) = Stau(Sixc(lp)) + (Sa(Sixc(lp))+(Sa(Six(lp))-Sa(Sixc(lp)))/2) *(abs(Sx(Sixc(lp))-Sx(Six(lp))));
    f = (numel(Six)+1 - lp)/numel(Six);
    waitbar(f,h);
end

Stau = Stau.*K^-1;

close(h);
end