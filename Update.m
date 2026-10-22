function [p, g] = Update(p,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE elevations using operator splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Order of the operators. With the implicit
% fluvial solvers (the default, and required for sediment) the time step is
% p.dtmax, known in advance, so the river step comes last: uplift, then
% hillslope transport, then rivers. The state at the end of each step (what is
% plotted and saved) is then the one right after the river step, so it does not
% show soil that hillslopes delivered to channels during the step and that the
% river has not yet had a chance to carry away, and no delivered soil is left
% pending for the next step. The hillslope step uses the channel mask from the
% previous river step (from TadpoleInitialize on the first step).
% The explicit fluvial solver sets the time step from the Courant number of
% the surface it erodes, so with it the river step comes first.

if p.doStreamPower && strcmp(p.fluvialSolver, 'explicit')

    % CHANNEL INCISION, forward-time explicit. StreamPower also
    % 1. adjusts the time step according to the stability factor
    % 2. records locations of channels
    [p, g] = StreamPower(p,g);

    % HILLSLOPE PROCESSES
    if p.doNLDiffusion
        g = NLDiff(p,g); % Nonlinear hillslope diffusion
    end

    % SOURCE TERMS / PERTURBATIONS
    g = Source(p,g);

else

    % SOURCE TERMS / PERTURBATIONS
    g = Source(p,g);

    % HILLSLOPE PROCESSES
    if p.doNLDiffusion
        g = NLDiff(p,g); % Nonlinear hillslope diffusion
    end

    % CHANNEL INCISION, and SEDIMENT ENTRAINMENT AND DEPOSITION if p.doSediment
    % (the fluvial solvers replace g.Budget; keep this step's hillslope terms)
    if p.doStreamPower
        hillBudget = struct();
        if isfield(g,'Budget')
            for f = fieldnames(g.Budget)'
                if startsWith(f{1}, 'hill'), hillBudget.(f{1}) = g.Budget.(f{1}); end
            end
        end
        [p, g] = FluvialImplicit(p,g); % Backward-time implicit, with sediment if p.doSediment (C MEX solver mexSPACE)
        for f = fieldnames(hillBudget)'
            g.Budget.(f{1}) = hillBudget.(f{1});
        end
    end

end
