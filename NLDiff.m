function g = NLDiff(p,g)

% Hillslope transport over one time step p.dt: nonlinear diffusion with a
% critical slope Sc (Roering et al., 1999), with the diffusivity approximated
% by the series of Ganti et al. (2012), D = Knl*(1 + s + ... + s^N), s =
% (S/Sc)^2, N = p.hillN; linear diffusion if p.Sc = Inf.

Zstart = g.Z;
Deliv = zeros(p.K,p.J); % thickness of material delivered to cells that don't diffuse (fixed points, and channels if p.doChannelDiffusion = 0)

% Solver for this step. With p.hillSolver = 'auto' (the default): linear
% diffusion (Sc = Inf) by the implicit linear solver, otherwise the explicit
% solver when Phi = dt*max(Knl)/min(dx,dy)^2 is at most p.hillPhiSwitch
% (where it is cheaper), and the implicit solver (Picard iteration) above it; both solve
% the same law.
solver = p.hillSolver;
if strcmp(solver,'auto')
    Phi = p.dt*max(p.Knl(:))/min(p.dx,p.dy)^2;
    if isinf(p.Sc)
        solver = 'linear';
    elseif Phi > p.hillPhiSwitch
        solver = 'implicit';
    else
        solver = 'explicit';
    end
end
g.hillSolverUsed = solver;

if strcmp(solver,'implicit')

    % Fully implicit (backward Euler), solved by Picard iteration with
    % Anderson acceleration in C. Same law as the explicit solver; stable
    % for any time step
    [g.Z, Deliv, g.hillInfo] = HillslopePicard(p,g);
    g.hillRate = g.hillInfo.rate; % the starting guess for the next step

elseif strcmp(solver,'linear')

    % Implicit linear diffusion (p.Sc = Inf), in C
    [g.Z, Deliv, g.hillInfo] = HillslopeLinear(p,g);
    g.hillRate = g.hillInfo.rate; % the starting guess for the next step

else % explicit, in sub-steps for stability

    % Algorithm for stable sub-timesteps:

    % Get the time step requested by the user and calculate the amount of that step remaining.

    dttodo = p.dt;

    % While there is still time remaining to do,
    while dttodo > 0

        % Get the RHS and dtstab. Make sure dtstab is positive, throw an error if it is not
        [RHS, dtstab, RHSsink] = AndrewsHanks(p,g);
        if dtstab <= 0
            error(['The calculated maximum stable time step of ' num2str(dtstab) ' for the nonlinear diffusion term is invalid.']);
        end

        % take a sub-step of p.hillStabilityFactor (default 0.8) times the maximum
        % stable time step dtstab, or the time remaining if that is shorter
        if p.hillStabilityFactor*dtstab < dttodo
            dtsub = p.hillStabilityFactor*dtstab;
        else
            % take a step of dtsub = the time remaining
            dtsub = dttodo;
        end

        g.Z = g.Z + dtsub*RHS; % take the step
        Deliv = Deliv + dtsub*RHSsink;
        dttodo = dttodo - dtsub; % update time remaining to do
    end

end

% Material delivered to cells that don't diffuse. Channel cells keep a fixed
% elevation during the hillslope step; with sediment, what they receive is
% stored in g.SoilDeliv (solid volume per cell) and enters the river's
% sediment flux as a lateral input in the fluvial step that follows, where
% SPACE deposition decides how much settles and how much travels on.
% Without sediment, it is removed, as if carried away by the river. Material
% delivered to fixed points leaves the model. A channel cell that stands
% above a hillslope neighbour passes no hillslope flux to it (closed face),
% so channel cells only receive material.
if p.doChannelDiffusion
    chanSink = false(p.K,p.J);
else
    chanSink = g.C==1 & g.Channels;
end
fixedSink = g.C == 0;

% fraction of each cell inside the domain when reflected across mirror
% boundaries, for the budget
omega = ones(p.K,p.J);
if p.bvec(1) == 1, omega(:,1)   = omega(:,1)/2;   end
if p.bvec(2) == 1, omega(:,end) = omega(:,end)/2; end
if p.bvec(3) == 1, omega(1,:)   = omega(1,:)/2;   end
if p.bvec(4) == 1, omega(end,:) = omega(end,:)/2; end
wa = omega*p.dx*p.dy;

if p.doSediment
    % Hillslopes are soil-mantled: material moved by diffusion is soil with
    % the same porosity as river sediment. A gain adds to the sediment layer;
    % a loss removes sediment first, then bedrock (soil production).
    % Soil delivered to channel sinks: solid volume per cell (full cell, like
    % the fluxes Qs), passed to the fluvial step (FluvialImplicit).
    toRiver = chanSink & Deliv > 0;
    g.SoilDeliv = zeros(p.K,p.J);
    g.SoilDeliv(toRiver) = (1-p.phi)*p.dx*p.dy*Deliv(toRiver);
    gain = g.Z - Zstart;
    dH = gain;
    loss = gain < 0;
    dH(loss) = -min(g.H(loss), -gain(loss));
    dR = gain - dH;
    g.H = g.H + dH;
    g.R = g.R + dR;
    % budget of this diffusion step (solid volumes)
    g.Budget.hillProduced = (1-p.phi)*sum(wa(:).*(-dR(:)));           % soil made from bedrock
    g.Budget.hillStorage = (1-p.phi)*sum(wa(:).*dH(:));               % change in stored sediment and soil
    g.Budget.hillExported = (1-p.phi)*sum(wa(fixedSink).*Deliv(fixedSink)); % delivered to fixed points
    g.Budget.hillToChannels = (1-p.phi)*sum(wa(toRiver).*Deliv(toRiver)); % delivered to channel sinks (the fluvial step's lateral input)
    g.Budget.hillResidual = g.Budget.hillProduced - g.Budget.hillStorage - g.Budget.hillExported ...
        - g.Budget.hillToChannels;
else
    % bulk volume removed at channels and fixed points
    g.Budget.hillExported = sum(wa(fixedSink | chanSink).*Deliv(fixedSink | chanSink));
    g.Budget.hillToChannels = sum(wa(chanSink & Deliv > 0).*Deliv(chanSink & Deliv > 0)); % delivered to channel sinks
end
