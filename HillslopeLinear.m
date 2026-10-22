function [Z, Deliv, info] = HillslopeLinear(p,g)

% Linear hillslope diffusion (p.Sc = Inf), q = Knl S, over one time step
% p.dt, implicitly by mexHillLinear: Peaceman-Rachford ADI when channels
% diffuse (p.doChannelDiffusion = 1) and Phi = dt*max(Knl)/min(dx,dy)^2 <= 0.2,
% and backward Euler (conjugate gradients) otherwise. Next to cells that
% don't diffuse (channels when they are sinks, fixed points and fixed edges),
% ADI drains too much soil into them every step, so backward Euler is used
% whenever channels are sinks. Stable for any time step; volume is conserved
% to round-off.
%
% Outputs:
%   Z      new elevations
%   Deliv  thickness of soil delivered to cells that don't diffuse (zero
%          elsewhere), as in NLDiff
%   info   linearSolver ('ADI' or 'BE'), cgIterations, rate (rate of change,
%          the starting guess for the next step)

dt = p.dt;
if p.doChannelDiffusion
    dodiff = double(g.C == 1);
else
    dodiff = double(g.C == 1 & ~g.Channels);
end

Phi = dt*max(p.Knl(:))/min(p.dx,p.dy)^2;
useADI = p.doChannelDiffusion && Phi <= 0.2;
% Channel cells that don't diffuse receive soil from higher hillslope neighbours
% but give none to lower ones: those faces are closed (dodiff = 2)
if ~p.doChannelDiffusion
    dodiff(g.C == 1 & g.Channels) = 2;
end
Zguess = [];
if ~useADI && isfield(g,'hillRate') && isequal(size(g.hillRate), size(g.Z))
    Zguess = g.Z + dt*g.hillRate;
end
[Z, Deliv, info1] = mexHillLinear(g.Z, Zguess, dodiff, p.Knl, [p.dx p.dy dt 1e-4 1000 useADI], p.bvec);
info.linearSolver = 'BE';
if useADI, info.linearSolver = 'ADI'; end
info.cgIterations = info1(1);
info.rate = (Z - g.Z)/dt;
