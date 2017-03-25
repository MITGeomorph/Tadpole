function [p g] = Update(p,g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE elevations using operator splitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANNEL INCISION
if p.doStreamPower
    [p g] = FT(p,g); % Forward-time explicit
    % Note that FT performs two additional steps:
    % 1. adjusts time step according to the Courant number
    % 2. records locations of channels
end


% SOURCE TERMS / PERTURBATIONS 
g = Source(p,g);


% LANDSLIDING
if p.doLandslides
    g = Landslide(p,g);
end


% DIFFUSION
if p.doDiffusion
    if (p.doStreamPower && ~p.doChannelDiffusion) || (p.doStreamPower && p.doAdaptiveTimeStep)
        g = SetUpADI(p,g); % update ADI matrices to exclude channel points
        % in the future, we could do this more efficiently by constructing
        % archetypal ADI matrices at the beginning, and then zeroing out
        % rows or multiplying by dtnew/dtold as necessary, depending on 
        % where the channels are and what the new time step is
    end
    g = ADI(p,g); % alternating-direction implicit (Crank-Nicolson in 2D)
end