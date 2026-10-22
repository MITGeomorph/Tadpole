function [Z, Deliv, info] = HillslopePicard(p,g)

% Nonlinear hillslope diffusion over one time step p.dt, fully implicit
% (backward Euler), by Picard iteration with Anderson acceleration in C
% (mexHillPicard). The law is the explicit solver's (AndrewsHanks): D =
% Knl*(1 + s + ... + s^N), s = (S/Sc)^2, capped at its value at Sc,
% (N + 1)*max(Knl), with N = p.hillN (default 5), so the two solvers differ
% only in time stepping. Stable for any time step; volume is conserved to
% round-off. Channels may diffuse or be sinks.
%
% Outputs: Z (new elevations), Deliv (thickness delivered to cells that
% don't diffuse); info: iterations, converged, cgIterations, rate (the rate
% of change, the starting guess for the next step).

if p.doChannelDiffusion
    dodiff = double(g.C == 1);
else
    dodiff = double(g.C == 1 & ~g.Channels);
end
% Channel cells that don't diffuse receive soil from higher hillslope neighbours
% but give none to lower ones: those faces are closed (dodiff = 2)
if ~p.doChannelDiffusion
    dodiff(g.C == 1 & g.Channels) = 2;
end
cap = p.hillN + 1; % D at Sc, as in AndrewsHanks

% starting guess: the previous step's rate of change
Zguess = [];
if isfield(g,'hillRate') && isequal(size(g.hillRate), size(g.Z))
    Zguess = g.Z + p.dt*g.hillRate.*(dodiff ~= 0);
end

% tolerances: Picard stops when the largest change between iterations is
% below 1e-4 of the largest change in the step; each linear solve to 1e-5
[Z, Deliv, i1] = mexHillPicard(g.Z, Zguess, dodiff, p.Knl, [p.dx p.dy p.dt p.Sc p.hillN cap 1e-4 100 6 1e-5], p.bvec);
info.iterations = i1(1);
info.converged = i1(2);
info.cgIterations = i1(3);
info.rate = (Z - g.Z)/p.dt;
