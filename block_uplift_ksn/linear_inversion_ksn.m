function [Ksn_Inv] = linear_inversion_ksn(DEM,varargin)
% X-z inversion for ksn after Goren et al. (2014) and Gallen (2018).
% Inputs:
%       (1) DEM clipped to a watershed as a TopoToolbox GRIDobj
%       (2) crita - critical drainage area for channel head initiation in
%           m^2 (default = 1e6)
%       (3) mn - the m over n ratio (concavity) of river system (default =
%           0.5)
%       (4) chi_inc - chi binning increment size solved for in inversion
%           (default = 1)
%       (5) gam - Gamma is a dampening coefficient - higher values result 
%           in more smoothing. If a value of -1 is used, this code will use
%           an L-curve optimazation scheme to find the proper coefficient.
%           (default = 10)
%       (6) flowOption: string of 'fill' or 'carve' for flow routing.
%          (optional) {default --> empty assumes that DEM is already 
%          filled or carved}
%
% Outputs:
%       A single matlab data structure Ksn_Inv that contains
%       (1) Ksn_Inv.A - forward model matrix
%       (2) Ksn_Inv.ksn - chi-average ksn
%       (3) Ksn_Inv.S - topotoolbox STREAMobj
%       (4) Ksn_Inv.Sz - elevation of network relative to outlet elevation
%       (5) Ksn_Inv.Schi - chi of the stream network
%       (6) Ksn_Inv.chi_steps - increments of chi used in inversion
%       (7) Ksn_Inv.res - resolution matrix
%       (8) Ksn_Inv.gam - gamma 
%
% Author: Sean F. Gallen
% Date Modified: 09/01/2023
% email: sean.gallen[at]colostate.edu

%%Parse Inputs
p = inputParser;         
p.FunctionName = 'linear_inversion_ksn';

% required inputs
addRequired(p,'DEM', @(x) isa(x,'GRIDobj'));

% optional inputs
addOptional(p,'crita', 1e6, @(x) isscalar(x));
addOptional(p,'mn', 0.5, @(x) isscalar(x));
addOptional(p,'chi_inc', 1, @(x) isscalar(x));
addOptional(p,'gam', 10, @(x) isscalar(x));
addOptional(p,'flowOption', []);

parse(p,DEM, varargin{:});
DEM   = p.Results.DEM;
crita = p.Results.crita;
mn = p.Results.mn;
chi_inc= p.Results.chi_inc;
gam = p.Results.gam;

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

Sz = DEM.Z(S.IXgrid)-min(DEM.Z(S.IXgrid));

% calculate chi
Sa = DA.Z(S.IXgrid);
a = Sa.^-mn;
Schi = cumtrapz(S,a);

% define desired time increments.
min_chi = floor(min(Schi));
max_chi = ceil(max(Schi));

if chi_inc > max_chi
    error(['Error: Your Tau increment is larger than the maximum Tau. ',...
        'The maximum chi is: ', num2str(max_chi),' years. Reduce the the increment',...
        'so that you have ~5 to 10 chi increments between 0 and the maximum chi value.']);
end

% make time vector
chi_steps = min_chi:chi_inc:max_chi; % chi steps
q = length(chi_steps);
n = length(Sz);  % number of nodes

% load the time matrix with data
A = zeros(n,q);
for i = 2:q
    inc_chi_t=chi_steps(i)-chi_steps(i-1);
    A(Schi >= chi_steps(i),i-1) = inc_chi_t;
end
for i = 1:n
    loc_sum = sum(A(i,:));
    loc = find(A(i,:)==0,1);
    A(i,loc) = Schi(i) - loc_sum; %this is the reminder
end

 % calculate ksn_pri
ksn_pri = ones(q,1)*(max(Sz)./max(Schi));

% make q by q identity matrix
I = eye(q,q);

if gam == -1
    % impose dampening (e.g. smoothing) parameter, (Gamma)
    Gam = logspace(-2,5,1000);
    
    % find the best-fit Gamma
    MisFit = nan(1,length(Gam));
    
    for i = 1:length(Gam)
        % now we can invert to find U following Goren from Tarantola, 1987
        ksn_mod = ksn_pri + (A'*A + Gam(i)^2*I)\A'*(Sz-A*ksn_pri);
        MisFit(i) = (1/(n-q))*sqrt(sum((Sz - A*ksn_mod).^2));
    end
    
    [ind,~] = turingPointFinder(1./Gam,MisFit);

    % define best-fit Gamma
    Gam = Gam(ind);
    disp(['Best fit Gamma = ',str2num(Gam)]);
else
    Gam = gam;
end

%% invert for uplift history with best damping factor
ksn_mod = ksn_pri + (A'*A + Gam^2*I)\A'*(Sz-A*ksn_pri);

% calculate the resolution
denom = (A'*A + Gam^2*I);
res = (denom\A')*A;

% Pack relevant results into structured array
Ksn_Inv.A = A;
Ksn_Inv.ksn_mod = ksn_mod;
Ksn_Inv.S = S;
Ksn_Inv.Schi = Schi;
Ksn_Inv.Sz = Sz;
Ksn_Inv.chi_steps = chi_steps;
Ksn_Inv.res = res;
Ksn_Inv.gam = gam;

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