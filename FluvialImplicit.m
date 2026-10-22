function [p, g] = FluvialImplicit(p,g)

% Implicit (backward Euler) solution of river incision and, if p.doSediment,
% sediment entrainment and deposition, for one time step of length p.dt.
%
% Physics (SPACE; Shobe et al., 2017), with Omega = (P*A)^m S^n:
%   bedrock erosion       Er = Kr Omega exp(-H/H*)        (bed lowering rate)
%   sediment entrainment  Es = Ks Omega (1 - exp(-H/H*))  (bed lowering rate)
%   deposition            D  = Vs Qs/Q                    (solid volume per unit bed area and time)
%   dR/dt = -Er,  dH/dt = D/(1-phi) - Es                  (uplift is applied in Source)
% Sediment fluxes Qs are solid volumes. A fraction Ff of eroded bedrock becomes
% wash load and leaves the model.
%
% Numerics: flow routing weights are
% computed once on the current surface and held fixed. Each cell's new
% elevation is found implicitly from its receivers' new elevations, visiting
% cells receivers first; the deposition term is coupled to the sediment flux
% from upstream by Gauss-Seidel iteration (Yuan et al., 2019). The solver is
% the C MEX file mexSPACE.

% flow routing on the current surface; weights are held fixed during the step
g = DrainageArea(p,g);

% cells where rivers erode (same rules as SPEqn). With p.channelThreshold =
% 'smooth' (default) every cell follows the erosion law with the smoothed
% thresholds p.Ecr and p.Ecs; with 'area' only cells with A > p.Ac erode.
[Ecr, Ecs] = ErosionThresholds(p);
if strcmp(p.channelThreshold, 'area')
    active = g.A > p.Ac;
else
    active = true(size(g.Z));
end
if p.flood
    erodes = g.C==1 & active & ~g.Flooded;
else
    erodes = g.C==1 & active & ~g.Minima;
end

if p.doSediment
    % Gauss-Seidel tolerance (Yuan et al., 2019)
    if isfield(p,'sedTol')
        tol = p.sedTol;
    else
        tol = 1e-3*max(abs(p.U(:)))*p.dt;
    end
    % Relaxation factor w for the inflow used for deposition. Where deposition is
    % strong (large Vs/P) and the time step long, a cell's inflow and the
    % elevations of the cell and its donors can feed back on each other so that
    % the iteration oscillates; the solver then reduces w during the step. The
    % value reached is carried to the next step (g.sedRelax) and allowed to
    % recover by a factor 1/0.7 per step.
    relax = 1;
    if isfield(g,'sedRelax'), relax = min(1, g.sedRelax/0.7); end
    par = [p.dx p.dy p.dt p.m p.n p.phi p.Ff p.Hstar tol p.sedMaxIter p.pitLimiter 1 relax 1 Ecr Ecs];
    % Lateral input (m^3/yr, solid volume per cell): the soil that the hillslope
    % step delivered to channel cells that do not diffuse (p.doChannelDiffusion =
    % 0), stored by NLDiff in g.SoilDeliv as a solid volume per cell. It is spread
    % over this step and enters each cell's sediment inflow, so that deposition
    % decides how much of it settles there and how much travels on. g.SoilDeliv
    % is then reset so that the same soil cannot be delivered twice.
    QsLat = [];
    if isfield(g,'SoilDeliv') && any(g.SoilDeliv(:))
        QsLat = g.SoilDeliv/p.dt;
    end
    g.SoilDeliv = zeros(p.K,p.J);
    [Z, R, H, Qs, QsIn, info, lakeFrac] = mexSPACE(g.Z, g.R, g.H, g.Qs, g.W, g.Stack, g.Q, g.Zf, g.Lake, ...
        double(erodes), double(g.C), p.Kr, p.Ks, p.Vs, par, p.bvec, QsLat);
    g.sedRelax = info(11);
    if ~info(9)
        warning('Tadpole:sedimentNotConverged', 'Sediment solver did not converge in %d iterations (max change %g m, tolerance %g m).', ...
            p.sedMaxIter, info(10), tol);
    end
    g.Z = Z; g.R = R; g.H = H; g.Qs = Qs; g.QsIn = QsIn;
    % mass budget of this step (solid volumes)
    g.Budget = struct('iterations', info(1), 'nLimited', info(2), 'produced', info(3), 'washload', info(4), ...
        'storage', info(5), 'exported', info(6), 'deficit', info(7), 'residual', info(8), 'relax', info(11), ...
        'lateral', info(12), 'lakeFrac', lakeFrac);
    g.Budget.newtonIterations = info(13); % mean Newton iterations per cell solve with thresholds (diagnostic)
else
    % stream power only: one sweep, receivers first
    par = [p.dx p.dy p.dt p.m p.n 0 0 1 0 1 p.pitLimiter 0 1 0 Ecr 0];
    [Z, ~, ~, ~, ~, info] = mexSPACE(g.Z, [], [], [], g.W, g.Stack, g.Q, g.Zf, g.Lake, ...
        double(erodes), double(g.C), p.Kr, [], [], par, p.bvec);
    g.Z = Z;
    g.Budget = struct('iterations', 1, 'nLimited', info(2));
end

% channel mask for the next hillslope step (ChannelMask: with 'smooth', Omega >
% Omega_c on the new surface; with 'area', A > p.Ac)
g.Channels = ChannelMask(p, g);
