function [p g] = FT(p,g)

% calculate RHS of the NKWE, recording the Courant number and which points 
% exceed the channel incision threshold
[RHS g.Channels Courant] = NKWE(g,p);

% adjust the time step based on the specified maximum Courant number and
% maximum time step
if p.doAdaptiveTimeStep
    p.dt = min( [p.dtmax p.Courant/Courant*p.dt]);
end

g.U = g.U + p.dt*RHS;