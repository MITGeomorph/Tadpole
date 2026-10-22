function TadpoleMovie(elevation,p,t,frame_interval,basename,sediment)

% TadpoleMovie(elevation,p,t,frame_interval,basename,sediment)
%
% produces an image sequence from a saved Tadpole run: the saved elevation
% grids (elevation), parameters struct p and time vector t, drawing an image
% every frame_interval saved grids. The images are written to a new subfolder named
% basename (in the current folder), as basename0001.png, basename0002.png, ...
% The image sequence can then be used to produce an animation with, e.g.,
% Quicktime's "Load image sequence" function. Using the default Quicktime
% codec at 15 fps produces a smooth animation suitable for Powerpoint or
% Keynote.
%
% p.plottype = 'perspective' draws an oblique 3-D view (PerspectiveView) with one
% vertical exaggeration for all frames (p.plotVE, or chosen from the largest relief).
% For runs with sediment, p.plottype = 'sediment' draws shaded relief colored
% by sediment thickness; pass the saved sediment thicknesses as the 6th argument.

[p.K,p.J,p.nSaved] = size(elevation);

if ~isfield(p,'plottype')
    p.plottype = 'elevation';
end

if isempty(t)
    t = (0:p.nSaved-1)*p.saveint;
end

p.zmin = min(elevation(:));
p.zmax = max(elevation(:));
if strcmp(p.plottype,'perspective') && ~isfield(p,'plotVE')
    % one vertical exaggeration for all frames, so that the growth of relief
    % shows: the largest relief of any frame spans 1/8 of the domain
    rel = max(max(elevation,[],[1 2]) - min(elevation,[],[1 2]));
    p.plotVE = max(1, max(p.dx*(p.J-1), p.dy*(p.K-1))/8/max(rel, eps));
end
if strcmp(p.plottype,'sediment')
    if nargin < 6
        error('Tadpole:movie', 'p.plottype = ''sediment'' needs the saved sediment thicknesses (sediment) as the 6th argument.');
    end
    p.Hmin = 0;
    p.Hmax = max([sediment(:); eps]); % sediment thickness color scale, fixed for all frames
end

% the images go in a subfolder named basename
if ~exist(basename, 'dir')
    mkdir(basename);
end
[~, name] = fileparts(basename);
framefile = fullfile(basename, name);

% loop through the frames
framenum = 0;

for k = 1:frame_interval:p.nSaved

    framenum = framenum+1;
    
    g.Z = elevation(:,:,k);
    if strcmp(p.plottype,'sediment'), g.H = sediment(:,:,k); end
    p.t = t(k);
    
    if k==1
        [p,g] = SetUpPlotMovie(p,g);
    else
        DrawPlotMovie((k-1)*p.saveint,p,g);
    end
   
    
    % print(fig, '-dtiff', '-r0', [framefile num2str(framenum, '%04d')]); % tiff, screen resolution
    
    % hide the axes toolbars, which MATLAB warns about when printing a figure
    ax = findall(p.fighandle, 'Type', 'axes');
    for j = 1:numel(ax)
        if isprop(ax(j), 'Toolbar') && ~isempty(ax(j).Toolbar), ax(j).Toolbar.Visible = 'off'; end
    end
    print(p.fighandle, '-dpng', [framefile num2str(framenum, '%04d')]); % png
        
end
    
% close the figure
close(p.fighandle);
