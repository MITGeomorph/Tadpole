function TadpoleMovie(output,p,t,frame_interval,basename)

% TadpoleMovie(output,p,t,frame_interval,basename)
%
% produces an image sequence from a saved Tadpole run using saved elevation
% grids in output, parameters struct p, time vector t, drawing an image every
% frame_interval saved grids with image file name root specified by basename.
% The image sequence can then be used to produce an animation with, e.g.,
% Quicktime's "Load image sequence" function. Using the default Quicktime
% codec at 15 fps produces a smooth animation suitable for Powerpoint or
% Keynote.

[p.K,p.J,p.N] = size(output);

if ~isfield(p,'plottype')
    p.plottype = 'elevation';
end

if isempty(t)
    t = (0:p.N-1)*p.saveint;
end


% loop through the frames
framenum = 0;

for n=1:frame_interval:p.N

    framenum = framenum+1;
    
    g.U = output(:,:,n);
    p.t = t(n);
    
    if n==1
        [p,g] = SetUpPlot(p,g);
    else
        DrawPlot((n-1)*p.saveint,p,g);
    end
   
    
    % print(fig, '-dtiff', '-r0', [basename num2str(framenum, '%04d')]); % tiff, screen resolution
    
    print(p.fighandle, '-dpng', [basename num2str(framenum, '%04d')]); % png
        
end
    
% close the figure
close(p.fighandle);
