function [p,g] = TadpoleRun(p,g)

% TadpoleRun.m
%
% Performs main iteration loop of Tadpole

n=0;

while p.t < p.tf
    
    n = n + 1;  
    
    %%%%%%%%%%%%%%%%%%%%%% UPDATE ELEVATIONS %%%%%%%%%%%%%%%%%%%%%%%

    [p g] = Update(p,g);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%% INCREMENT TIME %%%%%%%%%%%%%%%%%%%%%%%%%
    
    p.t = p.t + p.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if it is time to redraw the plot, do so.
    if p.doDrawPlot
        if ~rem(n,p.plotint)
            DrawPlot(n,p,g)
        end        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if p.doSaveOutput
        if ~rem(n,p.saveint)
            p.lastsave = n/p.saveint + 1;
            g.output(:,:,p.lastsave) = g.U;
            g.t(p.lastsave) = p.t;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

p.iterations = n;