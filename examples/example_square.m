% example_square.m
%
% sample script for running Tadpole


clear

%% SET PARAMETERS %%

% --------------- space and time resolution ------------------------------- 

p.Nx = 200;                 %     p.Nx             Number of grid points in x direction
p.Ny = 200;                 %     p.Ny             Number of grid points in y direction
p.dx = 250;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 250;                 %     p.dy             Grid spacing in the y direction (m)

p.doAdaptiveTimeStep = 1;   % p.doAdaptiveTimeStep Turn adaptive time step based on Courant number on (1) or off (0). If set to off, time step is p.dtmax
p.dtmax = 10000;            %     p.dtmax          maximum time step (yr)
p.Courant = 0.9;            %     p.Courant        maximum Courant number
p.tf = 1e6;                 %     p.tf             Total time of the simulation (yr)


% ----- boundary conditions, source terms, and flow routing ---------------

p.bdy.left  = 'fixed';      %     p.bdy            a struct with fields 'left' 'right' 'lower' 'upper'
p.bdy.right = 'fixed';      %                      specifying boundary condition:
p.bdy.upper = 'fixed';      % 
p.bdy.lower = 'fixed';      %                      'fixed'    --> constant elevation (derivatives set equal to zero)
                            %                      'mirror'   --> "mirror" boundary (centered differencing using boundary
                            %                                     point and one interior point)
                            %                      'periodic' --> periodic boundary (centered differencing using boundary
                            %                                     point, one interior point, and one
                            %                                     point from opposite boundary)

p.E = 1e-3;                 %     p.E              Rate of surface uplift or base level lowering (m/yr)

p.routing = 'Dms';          %     p.routing        Choose which flow routing method to use: 'D8' (steepest descent), 'Dinf' (Tarboton's D-infinity), or 'Dms' (multi-slope)
p.flood = 1;                %     p.flood          1=route flow through local minima, 0=don't

p.F = zeros(p.Ny,p.Nx);     %     p.F              Optional matrix of fixed points, in addition to boundary conditions above. Points with p.F == 1 will have constant elevation

                            
% ------------------ plotting and output ----------------------------------

p.doDrawPlot = 1;           %     p.doDrawPlot     Display solution as the run progresses
p.plotint = 200;            %     p.plotint        Plot will be redrawn every plotint iterations
p.plottype = 'color shade'; %     p.plottype       Options are 'perspective', 'mesh', 'drainage area', 'curvature', 'elevation', 'contour', 'shade', 'color shade'
p.zfactor = 2;              %     p.zfactor        Vertical exaggeration (optional, only applies to some types of plots)
                            %
p.doSaveOutput = 0;         %     p.SaveOutput     Save model output to a .mat file
p.saveint = 100;            %     p.saveint        Elevation grid will be saved every saveint iterations
p.runname = 'example_square'; %     p.runname:       Character string naming the run. If specified 
                            %                      (and if p.saveint~=0), the parameters and elevations at each 
                            %                      save interal will be saved in a binary .MAT file called <runname>.mat
                           
 
                            
% ------------------ hillslope processes ----------------------------------                            
                                                    
p.doDiffusion = 0;          %     p.doDiffusion    Turn hillslope diffusion on (1) or off (0)
p.D = 0.005;                %     p.D              Hillslope diffusivity (m^2/yr)
                            %
p.doLandslides = 1;         %     p.doLandslides   Turn landslides on (1) or off (0)
p.Sc = 0.6;                 %     p.Sc             Critical slope (m/m)


% ---------------- bedrock channel incision -------------------------------                           

p.doStreamPower = 1;        %     p.doStreamPower  Turn bedrock channel incision on (1) or off (0)
p.doChannelDiffusion = 0;   %     p.doChannelDiffusion Turn diffusion in channels on (1) or off (0)
p.Kf = 7.5e-3;              %     p.Kf             Coefficient in stream power incision law (kg m^(1+2m) yr^-2)
p.m = 0.5;                  %     p.m              Drainage area exponent in stream power law
p.w = 1.0;                  %     p.w              Slope exponent in stream power law
p.Kw = 1;                   %     p.Kw             Coefficient relating channel width to drainage area
p.wexp = 0;                 %     p.wexp           Exponent relating channel width to drainage area
p.thetac = 0;               %     p.thetac         Threshold for fluvial incision



% ------------------ initial conditions -----------------------------------                           

p.beta = 1.5;                 %     p.beta           Negative slope of the power spectrum. 0 = white noise, more positive values are "redder" (more variance at longer wavelengths)
p.variance = 1;             %     p.variance       Variance of elevation (m^2)
p.periodic = 1;             %     p.periodic       Elevations will be periodic at the boundaries (1) or not (0, default)

%% CREATE INITIAL SURFACE %%

init = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);

% level the boundaries and set minimum elevation to zero
zmin = min(init(:));
init(:,1) = zmin; % left
init(:,end) = zmin; % right
init(1,:) = zmin; % upper
init(end,:) = zmin; % lower
init = init - zmin;


%% RUN THE MODEL %%

% run the model, storing the final elevation grid in solution
solution = Tadpole(init,p);
