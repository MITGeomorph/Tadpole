function [p, g] = TadpoleFinalize(p,g)

% TadpoleFinalize.m
%
% Performs finalization steps for Tadpole

if p.doSaveOutput % if we're saving results

    % add the final state, unless the last save interval already saved it
    if g.t(p.lastsave) < p.t
        g.savedZ(:,:,p.lastsave+1)=g.Z;
        if p.doSediment
            g.savedH(:,:,p.lastsave+1) = g.H;
        end
        g.t(p.lastsave+1)=p.t;
    end

    elevation = g.savedZ; % elevations at the saved times
    vars = {'p', 't', 'elevation'};
    if p.doSediment
        sediment = g.savedH; % sediment thickness at the same times
        vars = [vars {'sediment'}];
    end
    
    t = g.t;
    
    if isfield(p,'fighandle')
        fighandle = p.fighandle;
    
        p = rmfield(p,'fighandle'); % temporarily remove the figure handle to avoid a large output file
    
        save(p.runname, '-v7.3', vars{:});
    
        p.fighandle = fighandle; 
    else
        save(p.runname, '-v7.3', vars{:});
    end

    % take a snapshot of the final topography and save it as an image
    Snapshot(p,g);

end

if p.doDrawPlot % if we're plotting the solution
    DrawPlot(p.iterations,p,g)
end
