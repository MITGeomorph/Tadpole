function [p, g] = TadpoleInitialize(initial,p)

% TadpoleInitialize.m
%
% Performs initialization steps for Tadpole


% Tadpole 1 parameters that no longer exist (linear diffusion by ADI, landslides,
% the old incision threshold and the Courant number, now p.StabilityFactor)
old = {'doDiffusion','D','doLandslides','thetac','Courant'};
old = old(isfield(p, old));
if ~isempty(old)
    error('Tadpole:retiredParameter', ['This run script sets parameters of Tadpole 1 that Tadpole 2 does not use (p.%s). ' ...
        'See "Moving a run script from Tadpole 1" in README.md.'], strjoin(old, ', p.'));
end
% The channel width parameters are no longer used (the erosion law is E = Kr*(P*A)^m*S^n)
if isfield(p,'Kw') || isfield(p,'wexp')
    warning('Tadpole:retiredParameter', 'p.Kw and p.wexp are no longer used and are ignored.');
end

% Parameters renamed in Tadpole 2.0; catch run scripts that still use the old names
if isfield(p,'E') && ~isfield(p,'U')
    error('Tadpole:renamedParameter', 'The uplift rate is now called p.U (it used to be p.E). Please rename p.E to p.U in the run script.');
end
if isfield(p,'Kf') && ~isfield(p,'Kr')
    error('Tadpole:renamedParameter', 'The bedrock erodibility is now called p.Kr (it used to be p.Kf). Please rename p.Kf to p.Kr in the run script.');
end
% The slope exponent of the stream power law is now p.n (it used to be p.w)
if isfield(p,'w')
    if ~isfield(p,'n'), p.n = p.w; end
    warning('Tadpole:renamedParameter', 'The slope exponent is now called p.n (it used to be p.w). This run uses p.n = %g.', p.n);
end

%%%%%%%%%%%%%%%%%%
% TIME VARIABLES %
%%%%%%%%%%%%%%%%%%
                           
p.t = 0; % time in yr
p.dt = p.dtmax;
if ~isfield(p,'StabilityFactor')
    p.StabilityFactor = 0.5; % maximum Courant number of the explicit fluvial solver
end
if ~isfield(p,'hillStabilityFactor')
    p.hillStabilityFactor = 0.8; % sub-steps of the explicit hillslope solver: this fraction of the maximum stable time step
end
if ~(p.hillStabilityFactor > 0 && p.hillStabilityFactor <= 1)
    error('Tadpole:hillStabilityFactor', 'p.hillStabilityFactor must be in (0, 1]: the explicit hillslope solver is unstable above 1.');
end
% p.fixedTimeStep: 1 = every time step is p.dtmax; 0 (default) = the explicit
% fluvial solver may shorten a step to keep its Courant number below
% p.StabilityFactor. The implicit solvers always take steps of p.dtmax. The
% explicit hillslope solver always divides each step into stable sub-steps.
if isfield(p,'doAdaptiveTimeStep')
    if ~isfield(p,'fixedTimeStep'), p.fixedTimeStep = ~p.doAdaptiveTimeStep; end
    warning('Tadpole:renamedParameter', ['p.doAdaptiveTimeStep is now p.fixedTimeStep, with the opposite meaning: ' ...
        '1 = every time step is p.dtmax, 0 = the explicit fluvial solver may shorten steps for stability. ' ...
        'This run uses p.fixedTimeStep = %d.'], p.fixedTimeStep);
end
if ~isfield(p,'fixedTimeStep'), p.fixedTimeStep = 0; end
if isfield(p,'doHillSubsteps')
    warning('Tadpole:retiredParameter', 'p.doHillSubsteps is no longer used: the explicit hillslope solver always takes stable sub-steps.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL AND BOUNDARY CONDITIONS % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p.K, p.J] = size(initial);

g.Z = initial;

% Fluvial solver: 'implicit' (default; backward time, no stability limit,
% required for sediment) or 'explicit' (forward time, time step shortened to keep the Courant
% number below p.StabilityFactor unless p.fixedTimeStep = 1).
if ~isfield(p,'fluvialSolver')
    p.fluvialSolver = 'implicit';
end
if ~any(strcmp(p.fluvialSolver, {'explicit','implicit'}))
    error('Tadpole:fluvialSolver', 'p.fluvialSolver must be ''explicit'' or ''implicit''.');
end
% Hillslope solver: 'auto' (default: chosen each step by NLDiff; linear
% diffusion (Sc = Inf) by an implicit linear solver, otherwise 'explicit'
% when Phi = dt*max(Knl)/min(dx,dy)^2 <= p.hillPhiSwitch (default 0.3, where
% it is cheaper) and 'implicit' above), 'explicit' (forward time, sub-stepped
% for stability), or 'implicit' (backward Euler with the same law, by Picard
% iteration in C; stable for any time step, for long steps). p.hillN is the
% number of terms after the first in the series approximation of the
% nonlinear diffusivity (Ganti et al., 2012), default 5; the diffusivity is
% capped at its value at Sc, (N + 1)*max(Knl).
if ~isfield(p,'hillSolver')
    p.hillSolver = 'auto';
end
if strcmp(p.hillSolver,'picard'), p.hillSolver = 'implicit'; end % another name for the implicit solver
if ~any(strcmp(p.hillSolver, {'auto','explicit','implicit'}))
    error('Tadpole:hillSolver', 'p.hillSolver must be ''auto'', ''explicit'' or ''implicit''.');
end
if ~isfield(p,'hillPhiSwitch'), p.hillPhiSwitch = 0.3; end % crossover of the explicit and implicit costs
% Warn when the explicit solver chosen in the run script will be much slower
% than the default for these parameters (at the largest time step)
if p.doNLDiffusion
    PhiMax = p.dtmax*max(p.Knl(:))/min(p.dx,p.dy)^2;
    Nw = 5; if isfield(p,'hillN'), Nw = p.hillN; end % series terms (default set below)
    if strcmp(p.hillSolver,'explicit') && PhiMax > 2
        warning('Tadpole:slowHillSolver', ['The explicit hillslope solver with Phi = dt*Knl/dx^2 = %.3g may need up to ~%d sub-steps ' ...
            'per time step; p.hillSolver = ''auto'' (the default) would use the implicit solver for the same law, which is much faster here.'], ...
            PhiMax, ceil(4*(Nw + 1)/p.hillStabilityFactor*PhiMax));
    end
end
if ~isfield(p,'hillN')
    p.hillN = 5;
end
if p.hillN < 1 || p.hillN ~= round(p.hillN)
    error('Tadpole:hillN', 'p.hillN must be a positive integer.');
end

if ~isfield(p,'pitLimiter')
    p.pitLimiter = 1; % the implicit solver may not lower a cell below its lowest receiver
end

% Where rivers erode. p.channelThreshold = 'smooth' (default): every cell
% follows the erosion law with smoothed thresholds, F(K Omega; Ec) = K Omega -
% Ec (1 - exp(-K Omega/Ec)), with the threshold erosion rates p.Ecr (bedrock)
% and p.Ecs (sediment) in m/yr, and channels (for hillslope transport) are
% cells where the stream power exceeds its threshold. 'area': rivers erode
% only where A > p.Ac (the rule of Tadpole 1). A run script that sets p.Ac
% but not p.Ecr runs with 'area'.
if ~isfield(p,'channelThreshold')
    if isfield(p,'Ecr') || ~p.doStreamPower
        p.channelThreshold = 'smooth';
    elseif isfield(p,'Ac')
        p.channelThreshold = 'area';
        warning('Tadpole:areaThreshold', ['p.Ac is set but p.Ecr is not, so rivers erode only where A > p.Ac ' ...
            '(p.channelThreshold = ''area''). The default is now a smoothed erosion threshold: set p.Ecr (and p.Ecs ' ...
            'with sediment), threshold erosion rates in m/yr (see the example run scripts and the README).']);
    else
        p.channelThreshold = 'smooth';
    end
end
if ~any(strcmp(p.channelThreshold, {'smooth','area'}))
    error('Tadpole:channelThreshold', 'p.channelThreshold must be ''smooth'' or ''area''.');
end
if p.doStreamPower
    if strcmp(p.channelThreshold, 'smooth')
        need = {'Ecr'}; if p.doSediment, need{end+1} = 'Ecs'; end
        for f = need
            if ~isfield(p, f{1})
                error('Tadpole:channelThreshold', ['p.channelThreshold = ''smooth'' needs p.%s, a threshold erosion rate (m/yr). ' ...
                    'A rough choice: channel heads at drainage area Ah and slope Sh erode at about K*(P*Ah)^m*Sh^n, so use ' ...
                    'that with K = p.Kr for p.Ecr and K = p.Ks for p.Ecs. 0 means no threshold.'], f{1});
            end
            if ~isscalar(p.(f{1})) || p.(f{1}) < 0
                error('Tadpole:channelThreshold', 'p.%s must be a non-negative scalar (m/yr).', f{1});
            end
        end
    elseif ~isfield(p,'Ac')
        error('Tadpole:channelThreshold', 'p.channelThreshold = ''area'' needs p.Ac, the drainage area threshold (m^2).');
    end
end

if p.doSediment
    if strcmp(p.fluvialSolver,'explicit')
        error('Tadpole:sediment', 'Sediment requires p.fluvialSolver = ''implicit''.');
    end
    if ~p.flood
        error('Tadpole:sediment', 'Sediment requires p.flood = 1 (sediment fills lakes up to the level of the filled surface).');
    end
    if ~p.doStreamPower
        error('Tadpole:sediment', 'Sediment requires p.doStreamPower = 1.');
    end
    for f = {'Ks','Hstar','Vs'}
        if ~isfield(p,f{1})
            error('Tadpole:sediment', 'Sediment requires p.%s (see the example run scripts).', f{1});
        end
    end
    if ~isfield(p,'phi'), p.phi = 0; end % porosity of sediment and soil
    if ~isfield(p,'Ff'), p.Ff = 0; end   % fraction of eroded bedrock lost as wash load
    if ~isfield(p,'sedMaxIter'), p.sedMaxIter = 200; end
    if ~isfield(p,'sedTol') && max(abs(p.U(:))) == 0
        error('Tadpole:sediment', ['With no uplift, the default tolerance of the sediment solver (1e-3*max(abs(p.U(:)))*p.dtmax) is zero. ' ...
            'Set p.sedTol to about 1e-3 times the typical elevation change per time step (see the example run scripts).']);
    end
    % initial sediment thickness (p.H0, scalar or matrix, optional); the bedrock surface lies below it
    if isfield(p,'H0')
        g.H = p.H0.*ones(p.K,p.J);
    else
        g.H = zeros(p.K,p.J);
    end
    g.R = g.Z - g.H;
    g.QsIn = zeros(p.K,p.J); % sediment flux entering each cell
    g.Qs = zeros(p.K,p.J);   % sediment flux leaving each cell (routed through the new weights, the starting guess for the next step)
    g.SoilDeliv = zeros(p.K,p.J); % soil delivered by hillslopes to channel cells that do not diffuse (solid volume), the next fluvial step's lateral input
end

[p, g] = BoundaryMat(p,g);

% Periodic edges must come in pairs (left with right, upper with lower).
% Flow routing, the sediment solver and hillslope diffusion all assume it:
% with one periodic edge, cells on that edge are connected to the opposite
% edge but not the other way, so hillslope fluxes are not conserved and the
% implicit hillslope matrix is not symmetric.
if xor(p.bvec(1)==2, p.bvec(2)==2) || xor(p.bvec(3)==2, p.bvec(4)==2)
    warning('Tadpole:unpairedPeriodic', ['Periodic boundaries must be paired (left with right, upper with lower), ' ...
        'but p.bdy has left = %s, right = %s, upper = %s, lower = %s. Results will not be physically meaningful.'], ...
        p.bdy.left, p.bdy.right, p.bdy.upper, p.bdy.lower);
end

if isscalar(p.P)
    p.P = p.P*ones(p.K,p.J); % Flow routing functions expect P to be a matrix
end

% channel mask on the initial surface, used by the first hillslope step (later
% ones use the mask of the previous river step)
if p.doStreamPower && ~p.doChannelDiffusion
    g = DrainageArea(p,g);
    g.Channels = ChannelMask(p, g);
else
    g.Channels = zeros(p.K,p.J);
end

%%%%%%%%%%%%%%%%%
% SET UP OUTPUT % 
%%%%%%%%%%%%%%%%%

if p.doSaveOutput
    
    % assign a default name for the output if it was not specified
    if ~isfield(p,'runname')
        p.runname = 'default_run_name';
    end
    p.saveint = ceil(p.saveint);
    g.savedZ(:,:,1) = initial; % the initial surface
    if p.doSediment
        g.savedH(:,:,1) = g.H; % sediment thickness is saved too (bedrock is savedZ - savedH)
    end
    g.t = 0; % vector that will hold the times corresponding to the saved grids
    p.lastsave = 1; % index of the last saved grid
end


%%%%%%%%%%%%%%%%%%%%%
%   SET UP PLOT     % 
%%%%%%%%%%%%%%%%%%%%%

if p.doDrawPlot
    if strcmp(p.plottype,'sediment') && ~p.doSediment
        error('Tadpole:plottype', 'p.plottype = ''sediment'' requires p.doSediment = 1; use ''elevation'' or ''hillshade''.');
    end
    p.plotint = round(p.plotint);
    [p, g] = SetUpPlot(p,g);
end
