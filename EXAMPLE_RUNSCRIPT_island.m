% EXAMPLE_RUNSCRIPT_island.m
%
% Example run script for Tadpole: erosion of a volcanic island that is no
% longer growing.
%
% What it shows:
%   - fixed points inside the domain (p.F): cells at or below sea level keep
%     their elevation, so the coastline is the base level for the rivers
%   - spatially variable precipitation (p.P as a matrix): rainfall increases
%     toward the windward side, so valleys there grow faster
%   - no uplift (p.U = 0): the relief decays as rivers and hillslopes erode the
%     island
%   - a smoothed erosion threshold (p.Ecr), and channel cells that are sinks for
%     hillslope soil (p.doChannelDiffusion = 0)
%
% Run time: about 2 minutes on a laptop (without plotting).


clear

%% SET PARAMETERS %%

% --------------- space and time resolution -------------------------------

p.Nx = 400;                 %     p.Nx             Number of grid points in x direction
p.Ny = 400;                 %     p.Ny             Number of grid points in y direction
p.dx = 125;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 125;                 %     p.dy             Grid spacing in the y direction (m)

p.dtmax = 5000;             %     p.dtmax          Time step (yr). The implicit solvers (the default) are stable for any time step: choose it
                            %                      for accuracy. With the explicit fluvial solver, the maximum time step (see p.fixedTimeStep)
p.tf = 5e6;                 %     p.tf             Total time of the simulation (yr)
p.fixedTimeStep = 0;        %     p.fixedTimeStep  1 = every time step is p.dtmax. 0 = the explicit fluvial solver may shorten a time step to
                            %                      keep it stable. Matters only with p.fluvialSolver = 'explicit'


% ----- boundary conditions, source terms, and flow routing ---------------

p.bdy.left  = 'fixed';      %     p.bdy            A struct with fields 'left', 'right', 'upper' and 'lower' giving the boundary conditions:
p.bdy.right = 'fixed';      %                      'fixed'    --> constant elevation; water, sediment and soil leave the model here
p.bdy.upper = 'fixed';      %                      'mirror'   --> reflection across the edge: nothing crosses it (a drainage divide)
p.bdy.lower = 'fixed';      %                      'periodic' --> the edge is joined to the opposite edge (set both edges of the pair)

p.U = 0;                    %     p.U              Rate of rock uplift relative to the fixed boundaries (m/yr). Can be a scalar or a matrix of
                            %                      size (Ny,Nx)
% precipitation increasing linearly from 0.5 to 5 m/yr toward the upper right (the windward side)
minP = 0.5;
maxP = 5;
p.P = linearGradient(p.Ny, p.Nx, minP, maxP, 30); % p.P Mean annual precipitation rate (m/yr). Can be a scalar or a matrix of size (Ny,Nx)

p.routing = 'Dms';          %     p.routing        Flow routing method: 'D8' (steepest descent), 'Dinf' (Tarboton's D-infinity), or 'Dms'
                            %                      (multiple flow directions, in proportion to slope^p.mfdexp)
p.mfdexp = 2;               %     p.mfdexp         Slope exponent for 'Dms' routing (1 = in proportion to slope; larger values concentrate flow in
                            %                      the steepest direction). Ignored by 'D8' and 'Dinf'
p.flood = 1;                %     p.flood          1 = route flow through closed depressions (lakes, filled by Priority-Flood), 0 = don't.
                            %                      Must be 1 with sediment: sediment fills lakes up to their spill level

p.F = zeros(p.Ny,p.Nx);     %     p.F              Optional matrix of fixed points, in addition to the boundary conditions above. Points with
                            %                      p.F == 1 keep a constant elevation (they act like fixed boundaries)


% ------------------ plotting and output ----------------------------------

p.doDrawPlot = 1;           %     p.doDrawPlot     Display the solution as the run progresses
p.plotint = 20;             %     p.plotint        The plot is redrawn every p.plotint time steps
p.plottype = 'hillshade';   %     p.plottype       'elevation' (map colored by elevation), 'hillshade' (shaded relief colored by elevation),
                            %                      'sediment' (shaded relief colored by sediment thickness; requires p.doSediment = 1), or
                            %                      'perspective' (oblique 3-D view; p.plotVE optionally fixes the vertical exaggeration)
p.doSaveOutput = 0;         %     p.doSaveOutput   Save model output to <p.runname>.mat, and an image of the final surface to <p.runname>.png
p.saveint = 50;             %     p.saveint        The elevation grid is saved every p.saveint time steps
p.runname = 'example_island';  %     p.runname        Name of the run, used for the output files. The .mat file holds p, the times t, the
                            %                      elevations (elevation, Ny x Nx x number of saved grids) and, with sediment, the
                            %                      sediment thicknesses (sediment)


% ------------------ hillslope processes ----------------------------------

p.doNLDiffusion = 1;        %     p.doNLDiffusion  Turn hillslope transport on (1) or off (0)
p.Knl = 0.02;               %     p.Knl            Hillslope transport coefficient (m^2/yr). Can be a scalar or a matrix of size (Ny,Nx)
p.Sc = 1.0;                 %     p.Sc             Critical slope (m/m) of the nonlinear transport law. For linear diffusion, set p.Sc = Inf
% p.hillSolver = 'auto';    %     p.hillSolver     'auto' (default): each time step uses the cheaper of 'explicit' (forward time, sub-stepped
                            %                      for stability) and 'implicit' (backward Euler, stable for any time step), which solve the
                            %                      same law; with p.Sc = Inf, an implicit linear solver. 'explicit' or 'implicit' forces one
% p.hillN = 5;              %     p.hillN          Number of terms after the first in the series approximation of the nonlinear diffusivity
                            %                      (Ganti et al., 2012). Default 5


% ---------------- river incision -----------------------------------------

p.doStreamPower = 1;        %     p.doStreamPower  Turn river incision on (1) or off (0)
p.doChannelDiffusion = 0;   %     p.doChannelDiffusion Hillslope transport in channel cells on (1) or off (0). With 0, channel cells are sinks
                            %                      for hillslope soil; with sediment, that soil enters the river's sediment flux
p.Kr = 5e-7;                %     p.Kr             Bedrock erodibility in the stream power law E = Kr*(P*A)^m*S^n (m^(1-3m) yr^(m-1)). Can be
                            %                      a scalar or a matrix of size (Ny,Nx)
p.m = 0.5;                  %     p.m              Discharge exponent m in the stream power law
p.n = 1.0;                  %     p.n              Slope exponent n in the stream power law
p.Ecr = 1e-4;               %     p.Ecr            Threshold erosion rate for bedrock (m/yr). Rivers erode every cell at the rate F(Kr*Omega),
                            %                      with Omega = (P*A)^m*S^n and the smoothed threshold F(x) = x - Ecr*(1 - exp(-x/Ecr)) (Shobe et
                            %                      al., 2017); cells where Kr*Omega > Ecr count as channels. Rule of thumb: for channel heads at
                            %                      drainage area Ah and slope Sh, Ecr ~ Kr*(P*Ah)^m*Sh^n. 0 = no threshold (rivers erode
                            %                      everywhere). Thresholds matter most on fine grids: on coarse grids nearly every cell in
                            %                      steep terrain exceeds them
% p.channelThreshold = 'area'; p.Ac = 1e5; %       Alternative (the rule of Tadpole 1): rivers erode only where drainage area > p.Ac (m^2).
                            %                      Faster, much faster with sediment (with the smoothed thresholds every cell is an eroding
                            %                      cell, so a time step with sediment costs about 4-5 times as much), but hollows near
                            %                      channel heads whose drainage area falls below p.Ac fill with hillslope soil that rivers
                            %                      cannot remove
p.fluvialSolver = 'implicit'; %   p.fluvialSolver  'implicit' (default; backward time, stable for any time step; required for sediment) or
                            %                      'explicit' (forward time, time step limited by the Courant number)


% ------------------ sediment ---------------------------------------------

p.doSediment = 0;           %     p.doSediment     Turn sediment entrainment and deposition (SPACE; Shobe et al., 2017) on (1) or off (0).
                            %                      Requires p.flood = 1 and p.fluvialSolver = 'implicit'. The parameters below are used only
                            %                      with sediment
% p.H0 = 0;                 %     p.H0             Initial sediment thickness (m), scalar or matrix of size (Ny,Nx); default 0
p.Ks = 5e-7;                %     p.Ks             Sediment erodibility (same units as p.Kr). Can be a scalar or a matrix of size (Ny,Nx)
p.Ecs = 0;                  %     p.Ecs            Threshold erosion rate for sediment entrainment (m/yr), as p.Ecr but with p.Ks
p.Hstar = 1.0;              %     p.Hstar          Bed roughness length scale H* (m): where sediment is thinner than this, bedrock is partly exposed
p.Vs = 3.0;                 %     p.Vs             Effective settling velocity of sediment (m/yr), as in Davy and Lague (2009). Can be a scalar
                            %                      or a matrix of size (Ny,Nx). Larger values (relative to p.P) deposit more sediment and make
                            %                      the solver take more iterations per time step
p.phi = 0.3;                %     p.phi            Porosity of sediment and soil
p.Ff = 0;                   %     p.Ff             Fraction of eroded bedrock that becomes wash load and leaves the model
p.sedTol = 1e-4;            %     p.sedTol         Convergence tolerance of the sediment solver (m). Default 1e-3*max(abs(p.U(:)))*p.dtmax,
                            %                      0.1% of the largest uplift in one time step. With no uplift anywhere you must set it:
                            %                      about 1e-3 times the typical elevation change per time step

%% INITIAL CONDITIONS %%

% a volcanic shield: a cone with a slope of 0.07 and a radius of 80% of the half-width
xc = (p.Nx-1)/2+1; % center of the island
yc = (p.Ny-1)/2+1;
r = 0.8*(p.Nx-1)/2; % radius of the island (cells)
S = 0.07; % slope of the shield
[X,Y] = meshgrid(1:p.Nx,1:p.Ny);
D = sqrt((X-xc).^2 + (Y-yc).^2);
shield = (r - D)*p.dx*S;

% add random topography with a standard deviation of about 140 m
beta = 2.0;                 %     beta             Negative slope of the power spectrum. 0 = white noise, larger values are "redder"
variance = 20000;           %     variance         Variance of elevation (m^2)
periodic = 0;               %     periodic         Elevations periodic at the boundaries (1) or not (0)
rng('shuffle'); % a different random surface in every run; rng(1) instead gives the same surface every time
noise = RedNoise(p.Ny,p.Nx,beta,variance,periodic);
noise = noise - mean(noise(:));
init = shield + noise;

% hold the sea floor at sea level
sealevel = 0;
p.F(init <= sealevel) = 1;
init(p.F == 1) = sealevel;

%% RUN THE MODEL %%

% run the model, storing the final elevation grid in solution (or, if
% p.doSaveOutput = 1, all saved grids)
solution = Tadpole(init,p);

%% HELPER FUNCTIONS %%

function M = linearGradient(Ny, Nx, minval, maxval, az)
% Ny-by-Nx matrix graded linearly from minval to maxval in the direction az
% (degrees, counterclockwise from +x)

% coordinate grid centered on the image
[x, y] = meshgrid( (0:Nx-1) - (Nx-1)/2, (0:Ny-1) - (Ny-1)/2 );

% projection along direction az, normalized to [0,1]
proj = x * cosd(az) + y * sind(az);
pmin = min(proj(:));
pmax = max(proj(:));
if pmax == pmin
    t = zeros(size(proj));
else
    t = (proj - pmin) / (pmax - pmin);
end

% map to [minval, maxval]
M = minval + t * (maxval - minval);

end
