function [p g] = TadpoleInitialize(initial,p)

% TadpoleInitialize.m
%
% Performs initialization steps for Tadpole


%%%%%%%%%%%%%%%%%%
% TIME VARIABLES % 
%%%%%%%%%%%%%%%%%%
                           
p.t = 0; % time in yr
p.dt = p.dtmax;
if p.doAdaptiveTimeStep && ~isfield(p,'Courant')
    p.Courant = 0.5; % default maximum Courant number
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL AND BOUNDARY CONDITIONS % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[p.K p.J] = size(initial);

g.U = initial;
% g.Uprev = g.U; % elevations at previous time step

[p g] = BoundaryMat(p,g);

if p.doStreamPower && ~p.doChannelDiffusion
    g.Channels = ones(p.K,p.J);
    [p g] = MarkChannels(p,g);
else
    g.Channels = zeros(p.K,p.J);
end

%%%%%%%%%%%%%%%%%
% SET UP OUTPUT % 
%%%%%%%%%%%%%%%%%

if p.doSaveOutput
    
    % assign a default name for the output if it was not specified
    if ~isfield(p,'runname')
        p.runname = 'default_run_name';
    end
    p.saveint = ceil(p.saveint);
    g.output(:,:,1) = initial; % the initial surface
    g.t = 0; % vector that will hold the times corresponding to the saved grids
end


%%%%%%%%%%%%%%%%%%%%%
% MEMORY ALLOCATION % 
%%%%%%%%%%%%%%%%%%%%%
                           
% construct matrix operators used in the Alternating-Direction Implicit (ADI) scheme for linear diffusion
if p.doDiffusion
    g = SetUpADI(p,g);
end


%%%%%%%%%%%%%%%%%%%%%
%   SET UP PLOT     % 
%%%%%%%%%%%%%%%%%%%%%

if p.doDrawPlot
    p.plotint = round(p.plotint);
    [p g] = SetUpPlot(p,g);
end
