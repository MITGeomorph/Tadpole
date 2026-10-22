% EXAMPLE_RUNSCRIPT_sediment.m
%
% Example run script for Tadpole with sediment entrainment and deposition
% (SPACE; Shobe et al., 2017), after the foreland simulation of Yuan et al.
% (2019, their Figure 4).
%
% A mountain range 100 km wide uplifts on the left side of the domain. The
% foreland on the right, 50 km wide, does not uplift and drains to a fixed base
% level along the right edge. Rivers carry sediment eroded from the mountains
% into the foreland, where it builds alluvial fans that grow and merge into a
% sediment-mantled plain crossed by shifting channels. Plot type 'sediment'
% shows shaded relief colored by sediment thickness.
%
% What it shows:
%   - sediment (p.doSediment = 1): rivers erode bedrock and entrain and
%     deposit sediment
%   - spatially variable uplift (p.U as a matrix)
%   - multiple-flow routing ('Dms'), which lets channels on the fans divide
%   - a mirror boundary on the left edge and periodic upper and lower edges:
%     water and sediment leave only across the fixed right edge
%
% Run time: about 3.5 minutes on a laptop (without plotting).


clear

%% SET PARAMETERS %%

% --------------- space and time resolution -------------------------------

p.Nx = 300;                 %     p.Nx             Number of grid points in x direction
p.Ny = 200;                 %     p.Ny             Number of grid points in y direction
p.dx = 500;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 500;                 %     p.dy             Grid spacing in the y direction (m)

p.dtmax = 1e4;              %     p.dtmax          Time step (yr). The implicit solvers (the default) are stable for any time step: choose it
                            %                      for accuracy. With the explicit fluvial solver, the maximum time step (see p.fixedTimeStep)
p.tf = 6e6;                 %     p.tf             Total time of the simulation (yr)
p.fixedTimeStep = 0;        %     p.fixedTimeStep  1 = every time step is p.dtmax. 0 = the explicit fluvial solver may shorten a time step to
                            %                      keep it stable. Matters only with p.fluvialSolver = 'explicit'


% ----- boundary conditions, source terms, and flow routing ---------------

p.bdy.left  = 'mirror';     %     p.bdy            A struct with fields 'left', 'right', 'upper' and 'lower' giving the boundary conditions:
p.bdy.right = 'fixed';      %                      'fixed'    --> constant elevation; water, sediment and soil leave the model here
p.bdy.upper = 'periodic';   %                      'mirror'   --> reflection across the edge: nothing crosses it (a drainage divide)
p.bdy.lower = 'periodic';   %                      'periodic' --> the edge is joined to the opposite edge (set both edges of the pair)

% uplift only in the mountain block (the left 100 km, two thirds of the domain)
p.U = zeros(p.Ny,p.Nx);     %     p.U              Rate of rock uplift relative to the fixed boundaries (m/yr). Can be a scalar or a matrix of
                            %                      size (Ny,Nx)
p.U(:,1:round(2/3*p.Nx)) = 2e-4;
p.P = 1.0;                  %     p.P              Mean annual precipitation rate (m/yr). Can be a scalar or a matrix of size (Ny,Nx)

p.routing = 'Dms';          %     p.routing        Flow routing method: 'D8' (steepest descent), 'Dinf' (Tarboton's D-infinity), or 'Dms'
                            %                      (multiple flow directions, in proportion to slope^p.mfdexp)
p.mfdexp = 1;               %     p.mfdexp         Slope exponent for 'Dms' routing (1 = in proportion to slope; larger values concentrate flow in
                            %                      the steepest direction). Ignored by 'D8' and 'Dinf'
p.flood = 1;                %     p.flood          1 = route flow through closed depressions (lakes, filled by Priority-Flood), 0 = don't.
                            %                      Must be 1 with sediment: sediment fills lakes up to their spill level

p.F = zeros(p.Ny,p.Nx);     %     p.F              Optional matrix of fixed points, in addition to the boundary conditions above. Points with
                            %                      p.F == 1 keep a constant elevation (they act like fixed boundaries)


% ------------------ plotting and output ----------------------------------

p.doDrawPlot = 1;           %     p.doDrawPlot     Display the solution as the run progresses
p.plotint = 20;             %     p.plotint        The plot is redrawn every p.plotint time steps
p.plottype = 'sediment';    %     p.plottype       'elevation' (map colored by elevation), 'hillshade' (shaded relief colored by elevation),
                            %                      'sediment' (shaded relief colored by sediment thickness; requires p.doSediment = 1), or
                            %                      'perspective' (oblique 3-D view; p.plotVE optionally fixes the vertical exaggeration)
p.doSaveOutput = 0;         %     p.doSaveOutput   Save model output to <p.runname>.mat, and an image of the final surface to <p.runname>.png
p.saveint = 20;             %     p.saveint        The elevation grid is saved every p.saveint time steps
p.runname = 'example_sediment';  %     p.runname        Name of the run, used for the output files. The .mat file holds p, the times t, the
                            %                      elevations (elevation, Ny x Nx x number of saved grids) and, with sediment, the
                            %                      sediment thicknesses (sediment)


% ------------------ hillslope processes ----------------------------------

p.doNLDiffusion = 1;        %     p.doNLDiffusion  Turn hillslope transport on (1) or off (0)
p.Knl = 0.01;               %     p.Knl            Hillslope transport coefficient (m^2/yr). Can be a scalar or a matrix of size (Ny,Nx)
p.Sc = 1.0;                 %     p.Sc             Critical slope (m/m) of the nonlinear transport law. For linear diffusion, set p.Sc = Inf
% p.hillSolver = 'auto';    %     p.hillSolver     'auto' (default): each time step uses the cheaper of 'explicit' (forward time, sub-stepped
                            %                      for stability) and 'implicit' (backward Euler, stable for any time step), which solve the
                            %                      same law; with p.Sc = Inf, an implicit linear solver. 'explicit' or 'implicit' forces one
% p.hillN = 5;              %     p.hillN          Number of terms after the first in the series approximation of the nonlinear diffusivity
                            %                      (Ganti et al., 2012). Default 5


% ---------------- river incision -----------------------------------------

p.doStreamPower = 1;        %     p.doStreamPower  Turn river incision on (1) or off (0)
p.doChannelDiffusion = 1;   %     p.doChannelDiffusion Hillslope transport in channel cells on (1) or off (0). With 0, channel cells are sinks
                            %                      for hillslope soil; with sediment, that soil enters the river's sediment flux
p.Kr = 2e-5;                %     p.Kr             Bedrock erodibility in the stream power law E = Kr*(P*A)^m*S^n (m^(1-3m) yr^(m-1)). Can be
                            %                      a scalar or a matrix of size (Ny,Nx)
p.m = 0.4;                  %     p.m              Discharge exponent m in the stream power law
p.n = 1.0;                  %     p.n              Slope exponent n in the stream power law
p.Ecr = 0;                  %     p.Ecr            Threshold erosion rate for bedrock (m/yr). Rivers erode every cell at the rate F(Kr*Omega),
                            %                      with Omega = (P*A)^m*S^n and the smoothed threshold F(x) = x - Ecr*(1 - exp(-x/Ecr)) (Shobe et
                            %                      al., 2017); cells where Kr*Omega > Ecr count as channels. Rule of thumb: for channel heads at
                            %                      drainage area Ah and slope Sh, Ecr ~ Kr*(P*Ah)^m*Sh^n. 0 = no threshold (rivers erode
                            %                      everywhere). Thresholds matter most on fine grids: on coarse grids nearly every cell in
                            %                      steep terrain exceeds them.
% p.channelThreshold = 'area'; p.Ac = 1e5; %       Alternative (the rule of Tadpole 1): rivers erode only where drainage area > p.Ac (m^2).
                            %                      Faster, much faster with sediment (with the smoothed thresholds every cell is an eroding
                            %                      cell, so a time step with sediment costs about 4-5 times as much), but hollows near
                            %                      channel heads whose drainage area falls below p.Ac fill with hillslope soil that rivers
                            %                      cannot remove
p.fluvialSolver = 'implicit'; %   p.fluvialSolver  'implicit' (default; backward time, stable for any time step; required for sediment) or
                            %                      'explicit' (forward time, time step limited by the Courant number)


% ------------------ sediment ---------------------------------------------

p.doSediment = 1;           %     p.doSediment     Turn sediment entrainment and deposition (SPACE; Shobe et al., 2017) on (1) or off (0).
                            %                      Requires p.flood = 1 and p.fluvialSolver = 'implicit'. The parameters below are used only
                            %                      with sediment
% p.H0 = 0;                 %     p.H0             Initial sediment thickness (m), scalar or matrix of size (Ny,Nx); default 0
p.Ks = 2e-5;                %     p.Ks             Sediment erodibility (same units as p.Kr). Can be a scalar or a matrix of size (Ny,Nx)
p.Ecs = 0;                  %     p.Ecs            Threshold erosion rate for sediment entrainment (m/yr), as p.Ecr but with p.Ks
p.Hstar = 1.0;              %     p.Hstar          Bed roughness length scale H* (m): where sediment is thinner than this, bedrock is partly exposed
p.Vs = 1.0;                 %     p.Vs             Effective settling velocity of sediment (m/yr), as in Davy and Lague (2009). Can be a scalar
                            %                      or a matrix of size (Ny,Nx). Larger values (relative to p.P) deposit more sediment and make
                            %                      the solver take more iterations per time step. Here Vs = P, the
                            %                      case G = 1 of Yuan et al. (2019)
p.phi = 0;                  %     p.phi            Porosity of sediment and soil
p.Ff = 0;                   %     p.Ff             Fraction of eroded bedrock that becomes wash load and leaves the model
% p.sedTol = 1e-4;          %     p.sedTol         Convergence tolerance of the sediment solver (m). Default 1e-3*max(abs(p.U(:)))*p.dtmax,
                            %                      0.1% of the largest uplift in one time step. With no uplift anywhere you must set it:
                            %                      about 1e-3 times the typical elevation change per time step

%% INITIAL CONDITIONS %%

% a flat surface with random bumps up to 1 m high
rng('shuffle'); % a different random surface in every run; rng(1) instead gives the same surface every time
init = rand(p.Ny,p.Nx);
init(:,end) = 0; % base level

%% RUN THE MODEL %%

% run the model, storing the final elevation grid in solution (or, if
% p.doSaveOutput = 1, all saved grids)
solution = Tadpole(init,p);

%% MAKE FRAMES FOR A MOVIE %%

% (set p.doSaveOutput = 1 above first)
% basename = [p.runname '_movie']; % base name of the image files, and of the new subfolder they go in
% frame_interval = 2;          % draw every frame_interval-th saved grid
% load(p.runname)              % load the saved elevations (elevation), sediment thicknesses (sediment), parameters (p) and times (t)
% TadpoleMovie(elevation,p,t,frame_interval,basename,sediment) % write the frames as .png images
