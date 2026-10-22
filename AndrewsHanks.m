function [RHS, dtstab, RHSsink] = AndrewsHanks(p,g)

% calculates the RHS of the Nonlinear diffusion finite-difference approximation,
% using the method of Andrews & Hanks (1985), as described in Perron (2011).
% Records the maximum value of the nonlinear diffusivity on the grid.
% The face diffusivities come from HillslopeFaces.
%
% RHSsink is the rate at which cells that don't diffuse (fixed points, and
% channels if p.doChannelDiffusion = 0) receive material from their
% neighbors (the RHS before it is set to zero there), and zero elsewhere.

% Number of terms in the series approximation of Ganti et al. (2012)
if isfield(p,'hillN')
    N = p.hillN;
else
    N = 5;
end

% Cap on the nonlinear diffusivity, as a multiple of max(Knl): the value of the
% series at Sc, N + 1 (6 for N = 5), so above Sc transport continues with the
% diffusivity the law has at Sc and slopes far above Sc (left by the fluvial
% step) do not make the sub-timesteps very short. The implicit solver
% (HillslopePicard) uses the same cap.
capFactor = N + 1;

[Dnlr, Dnll, Dnlu, Dnld, dodiff, Dnlmax, z] = HillslopeFaces(p, g, g.Z, N, capFactor);

% Maximum stable time step of forward Euler for 2D diffusion on the five-point
% operator: dt <= 1/(2 D (1/dx^2 + 1/dy^2)), i.e. dx^2/(4 D) for square cells,
% half the 1D limit dx^2/(2 D). NLDiff takes sub-steps of a fraction
% p.hillStabilityFactor of it
dtstab = 1/(2*Dnlmax*(1/p.dx^2 + 1/p.dy^2));

% One-sided differences (z is padded according to the boundary conditions)
zxr = (z(2:end-1,3:end) - z(2:end-1,2:end-1)) / p.dx; % dz/dx, one-sided right difference (j+1 minus j)
zxl = (z(2:end-1,2:end-1) - z(2:end-1,1:end-2)) / p.dx; % dz/dx, one-sided left difference (j minus j-1)

zyu = (z(1:end-2,2:end-1) - z(2:end-1,2:end-1)) / p.dy; % dz/dy, one-sided up difference (i+1 minus i)
zyd = (z(2:end-1,2:end-1) - z(3:end,2:end-1)) / p.dy; % dz/dy, one-sided down difference (i minus i-1)

% Calculate the RHS
RHS = (Dnlr.*zxr - Dnll.*zxl) / p.dx + (Dnlu.*zyu - Dnld.*zyd) / p.dy;

% Record what the cells that don't diffuse receive, then make sure the RHS = zero there
RHSsink = RHS;
RHSsink(dodiff==1) = 0;
RHS(~dodiff) = 0;

% Channel cells that don't diffuse receive material from higher hillslope
% neighbours (RHSsink) but give none to lower ones: a face between a channel
% cell and a lower diffusing neighbour is closed.
% Mostly such neighbours are banks beside the next channel cell downstream of a
% channel that runs diagonally across the grid, so the flux would not be
% hillslope transport.
if ~p.doChannelDiffusion
    ch = g.C == 1 & g.Channels;
    inR =  Dnlr.*zxr/p.dx; inL = -Dnll.*zxl/p.dx; % flux into the cell across each face
    inU =  Dnlu.*zyu/p.dy; inD = -Dnld.*zyd/p.dy;
    RHSsink(ch) = max(inR(ch),0) + max(inL(ch),0) + max(inU(ch),0) + max(inD(ch),0);
    % diffusing cells: remove the inflow from higher channel neighbours
    chp = PadGrid(double(ch), p);
    cR = chp(2:end-1,3:end) == 1 & inR > 0; cL = chp(2:end-1,1:end-2) == 1 & inL > 0;
    cU = chp(1:end-2,2:end-1) == 1 & inU > 0; cD = chp(3:end,2:end-1) == 1 & inD > 0;
    dif = dodiff == 1;
    fix = dif & (cR | cL | cU | cD);
    RHS(fix) = RHS(fix) - inR(fix).*cR(fix) - inL(fix).*cL(fix) - inU(fix).*cU(fix) - inD(fix).*cD(fix);
end

end
