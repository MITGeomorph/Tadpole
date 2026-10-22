function [p, g] = StreamPower(p,g)

% calculate RHS of the NKWE, recording the Courant number and which points 
% exceed the channel incision threshold
[RHS, g] = SPEqn(g,p); % g includes g.Channels and g.Courant
% [RHS, g.Channels, Courant] = SPEqn(g,p);

% adjust the time step based on the specified maximum Courant number and
% maximum time step, unless the time step is fixed. With a fixed step, warn
% (once) if the Courant number exceeds 1: the explicit scheme is then unstable
if ~p.fixedTimeStep
    p.dt = min( [p.dtmax p.StabilityFactor/g.Courant*p.dt]);
elseif g.Courant > 1 && ~isfield(g,'courantWarned')
    warning('Tadpole:courant', ['With p.fixedTimeStep = 1, the Courant number of the explicit fluvial solver is %.3g (> 1), ' ...
        'so the solution may be unstable. Use a shorter p.dtmax, p.fixedTimeStep = 0, or p.fluvialSolver = ''implicit''.'], g.Courant);
    g.courantWarned = true;
end

g.Z = g.Z + p.dt*RHS;