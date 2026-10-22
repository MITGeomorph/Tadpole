function [Dr, Dl, Du, Dd, dodiff, Dnlmax, z] = HillslopeFaces(p, g, z, N, capFactor)

% Nonlinear hillslope diffusivities on the faces between cells, for the
% conservative finite-difference operator of Andrews & Hanks (1985), as
% described in Perron (2011), used by the explicit solver (AndrewsHanks):
%
%   RHS = (Dr.*zxr - Dl.*zxl)/dx + (Du.*zyu - Dd.*zyd)/dy
%
% Inputs:
%   z          elevations (K x J)
%   N          number of terms after the first in the series approximation
%              of Ganti et al. (2012): D = Knl*(1 + s + s^2 + ... + s^N),
%              s = (S/Sc)^2. D is positive and finite for all slopes.
%   capFactor  diffusivities are capped at capFactor*max(Knl) (Inf: no cap)
%
% Outputs:
%   Dr, Dl, Du, Dd  diffusivities on the right, left, up and down faces of each cell
%   dodiff          1 where cells diffuse, 0 at cells that don't (fixed
%                   points, and channels if p.doChannelDiffusion = 0)
%   Dnlmax          maximum cell-centre diffusivity
%   z               the elevations padded according to the boundary conditions

% Pad the grid according to boundary conditions p.bvec = [left right upper lower]
% (0 = fixed elevation, 1 = mirror, 2 = periodic)
z = PadGrid(z, p);

% Centered differences
zxc = (z(2:end-1,3:end) - z(2:end-1,1:end-2)) / (2*p.dx); % dz/dx, centered
zyc = (z(1:end-2,2:end-1) - z(3:end,2:end-1)) / (2*p.dy); % dz/dy, centered

% Calculate the nonlinear diffusivity at all points, using the Maclaurin
% series approximation of Ganti et al. (2012) with N terms after the first
SoverScsq = (zxc.*zxc + zyc.*zyc)/(p.Sc*p.Sc);
term = SoverScsq;
series = 1 + SoverScsq;
for k = 2:N
    term = term.*SoverScsq;
    series = series + term;
end
Dnl = p.Knl.*series;

% Set the nonlinear diffusivity to zero at all fixed points

% we have a matrix g.C that is 1 where uplift and erosion should occur,
% 0 where it shouldn't

if p.doChannelDiffusion
    dodiff = 1*g.C;
else
    dodiff = 1*(g.C & ~g.Channels);
end
Dnl(~dodiff) = 0;

% Optionally, set a max value for the nonlinear diffusivity (the explicit
% solver uses this because steep slopes can lead to very short sub-timesteps)
if isfinite(capFactor)
    Dnlcap = capFactor*max(p.Knl(:));
    Dnl(Dnl>Dnlcap) = Dnlcap;
end

Dnlmax = max(Dnl(:));

% Pad the nonlinear diffusivity matrix according to boundary conditions
Dnl = PadGrid(Dnl, p);

% Calculate the nonlinear diffusivity at half-points by averaging

Dr = 0.5*(Dnl(2:end-1,3:end) + Dnl(2:end-1,2:end-1)); % right (j+1 and j)
Dl = 0.5*(Dnl(2:end-1,2:end-1) + Dnl(2:end-1,1:end-2)); % left (j and j-1)

Du = 0.5*(Dnl(1:end-2,2:end-1) + Dnl(2:end-1,2:end-1)); % up (i+1 and i)
Dd = 0.5*(Dnl(2:end-1,2:end-1) + Dnl(3:end,2:end-1)); % down (i and i-1)

% Cells that don't diffuse (fixed points, and channels if diffusion is
% turned off in channels) act as sinks for hillslope sediment. Their
% diffusivity is zero, so averaging would halve the flux delivered across a
% face between a diffusing cell and a sink. On those faces, use the
% diffusivity of the diffusing cell instead (which equals the sum of the
% two diffusivities, since the sink's is zero).
chansink = PadGrid(1*(~dodiff), p); % cells that don't diffuse (sinks)
isdiff = PadGrid(dodiff, p); % cells that diffuse

oneside = (isdiff(2:end-1,3:end) & chansink(2:end-1,2:end-1)) | (chansink(2:end-1,3:end) & isdiff(2:end-1,2:end-1));
Dr(oneside) = 2*Dr(oneside); % right (j+1 and j)
oneside = (isdiff(2:end-1,2:end-1) & chansink(2:end-1,1:end-2)) | (chansink(2:end-1,2:end-1) & isdiff(2:end-1,1:end-2));
Dl(oneside) = 2*Dl(oneside); % left (j and j-1)
oneside = (isdiff(1:end-2,2:end-1) & chansink(2:end-1,2:end-1)) | (chansink(1:end-2,2:end-1) & isdiff(2:end-1,2:end-1));
Du(oneside) = 2*Du(oneside); % up (i+1 and i)
oneside = (isdiff(2:end-1,2:end-1) & chansink(3:end,2:end-1)) | (chansink(2:end-1,2:end-1) & isdiff(3:end,2:end-1));
Dd(oneside) = 2*Dd(oneside); % down (i and i-1)
