% example_Titan.m
%
% sample script for running Tadpole


clear

%% SET PARAMETERS %%

% --------------- space and time resolution ------------------------------- 

p.Nx = 200;                 %     p.Nx             Number of grid points in x direction
p.Ny = 200;                 %     p.Ny             Number of grid points in y direction
p.dx = 125;                 %     p.dx             Grid spacing in the x direction (m)
p.dy = 125;                 %     p.dy             Grid spacing in the y direction (m)

p.doAdaptiveTimeStep = 1;   % p.doAdaptiveTimeStep Turn adaptive time step based on Courant number on (1) or off (0). If set to off, time step is p.dtmax
p.dtmax = 1e6;              %     p.dtmax          maximum time step (yr)
p.Courant = 0.9;            %     p.Courant        maximum Courant number
p.tf = 1e6;                 %     p.tf             Total time of the simulation (yr)


% ----- boundary conditions, source terms, and flow routing ---------------

p.bdy.left  = 'periodic';   %     p.bdy            a struct with fields 'left' 'right' 'lower' 'upper'
p.bdy.right = 'periodic';   %                      specifying boundary condition:
p.bdy.upper = 'mirror';   % 
p.bdy.lower = 'fixed';   %                      'fixed'    --> constant elevation (derivatives set equal to zero)
                            %                      'mirror'   --> "mirror" boundary (centered differencing using boundary
                            %                                     point and one interior point)
                            %                      'periodic' --> periodic boundary (centered differencing using boundary
                            %                                     point, one interior point, and one
                            %                                     point from opposite boundary)

p.E = 1e-10;                %     p.E              Rate of surface uplift or base level lowering (m/yr)
p.routing = 'Dms';          %     p.routing        Choose which flow routing method to use: 'D8' (steepest descent), 'Dinf' (Tarboton's D-infinity), or 'Dms' (multi-slope)
p.flood = 1;                %     p.flood          1=route flow through local minima, 0=don't

p.F = zeros(p.Ny,p.Nx);     %     p.F              Optional matrix of fixed points, in addition to boundary conditions above. Points with p.F == 1 will have constant elevation
                            
% ------------------ plotting and output ----------------------------------

p.doDrawPlot = 1;           %     p.doDrawPlot     Display solution as the run progresses
p.plotint = 100;            %     p.plotint        Plot will be redrawn every plotint iterations
p.plottype = 'color shade';             %     p.plottype       1=perspective view, 2=drainage area map, 3=curvature map, 4=elevation map, 5=contour map, 6=shaded relief, 7=colored shaded relief
                            %
p.doSaveOutput = 0;         %     p.SaveOutput     Save model output to a .mat file
p.saveint = 100;              %     p.saveint        Elevation grid will be saved every saveint iterations
p.runname = 'example_Titan_linearblock';%     p.runname:       Character string naming the run. If specified 
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
p.Kf = 5e-6;                %     p.Kf             Coefficient in stream power incision law (kg m^(1+2m) yr^-2)
p.m = 0.5;                  %     p.m              Drainage area exponent in stream power law
p.w = 1.0;                  %     p.w              Slope exponent in stream power law
p.Kw = 125;                 %     p.Kw             Coefficient relating channel width to drainage area
p.wexp = 0;                 %     p.wexp           Exponent relating channel width to drainage area
p.thetac = 0;               %     p.thetac         Threshold for fluvial incision



% ------------------ initial conditions -----------------------------------                           

p.beta = 1.8;               %     p.beta           Negative slope of the power spectrum. 0 = white noise, more positive values are "redder" (more variance at longer wavelengths)
p.variance = 100;         %     p.variance       Variance of elevation (m^2)
p.periodic = 1;             %     p.periodic       Elevations will be periodic at the boundaries (1) or not (0, default)


%% CREATE INITIAL SURFACE %%

noise = RedNoise(p.Ny,p.Nx,p.beta,p.variance,p.periodic);

% % make lowest 10% of elevations fixed points
% z10 = prctile(noise(:),10);
% noise = noise - z10;
% noise(noise < 0) = 0; 
% p.F(noise <= 0) = 1;

% init = noise;

% add a background slope
bgslope = 0.05;
[~,Y] = meshgrid(p.dx*(0:p.Nx-1),p.dy*(p.Ny-1:-1:0));

init = noise + bgslope*Y;

% level lower boundary with elevation = 0
init(end,:) = 0;

%% RUN THE MODEL %%

% run the model, storing the final elevation grid in solution
solution = Tadpole(init,p);
