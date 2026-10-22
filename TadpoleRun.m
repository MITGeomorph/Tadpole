function [p,g] = TadpoleRun(p,g)

% TadpoleRun.m
%
% Performs main iteration loop of Tadpole

iter = 0;

while p.t < p.tf
    
    iter = iter + 1;  
    
    %%%%%%%%%%%%%%%%%%%%%% UPDATE ELEVATIONS %%%%%%%%%%%%%%%%%%%%%%%

    [p, g] = Update(p,g);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%% INCREMENT TIME %%%%%%%%%%%%%%%%%%%%%%%%%
    
    p.t = p.t + p.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % if it is time to redraw the plot, do so.
    if p.doDrawPlot
        if ~rem(iter,p.plotint)
            DrawPlot(iter,p,g)
        end        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if p.doSaveOutput
        if ~rem(iter,p.saveint)
            p.lastsave = iter/p.saveint + 1;
            g.savedZ(:,:,p.lastsave) = g.Z;
            if p.doSediment
                g.savedH(:,:,p.lastsave) = g.H;
            end
            g.t(p.lastsave) = p.t;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

p.iterations = iter;