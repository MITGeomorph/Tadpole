function TadpoleMovie(elevation,p,t,frame_interval,basename,varargin)

% TadpoleMovie(elevation,p,t,frame_interval,basename)
% TadpoleMovie(elevation,p,t,frame_interval,basename,sediment)
% TadpoleMovie(...,'fps',fps)
%
% produces an image sequence from a saved Tadpole run: the saved elevation
% grids (elevation), parameters struct p and time vector t, drawing an image
% every frame_interval saved grids. The images are written to a new subfolder named
% basename (in the current folder), as basename0001.png, basename0002.png, ...
%
% With 'fps', fps (frames per second, > 0), TadpoleMovie also makes a video
% of the images, basename.mp4 (next to the subfolder), as MPEG-4 (H.264),
% which plays in Powerpoint, Keynote, Quicktime and web browsers; 15 fps
% gives a smooth animation. Without it (or with fps = 0), only the images are
% written. MATLAB writes MPEG-4 only on Windows and macOS; elsewhere (Linux)
% the video is written as Motion JPEG AVI (basename.avi) instead, with a
% warning.
%
% p.plottype = 'perspective' draws an oblique 3-D view (PerspectiveView) with one
% vertical exaggeration for all frames (p.plotVE, or chosen from the largest relief).
% For runs with sediment, p.plottype = 'sediment' draws shaded relief colored
% by sediment thickness; pass the saved sediment thicknesses as the 6th argument.

% optional arguments: the sediment thicknesses, then name-value pairs
sediment = [];
if ~isempty(varargin) && isnumeric(varargin{1})
    sediment = varargin{1};
    varargin(1) = [];
end
fps = 0;
if mod(numel(varargin),2)
    error('Tadpole:movie', 'Optional arguments after the sediment thicknesses must be name-value pairs, such as ''fps'', 15.');
end
for j = 1:2:numel(varargin)
    switch lower(varargin{j})
        case 'fps'
            fps = varargin{j+1};
        otherwise
            error('Tadpole:movie', 'Unknown option ''%s''. The only option is ''fps''.', varargin{j});
    end
end
if ~isscalar(fps) || ~isnumeric(fps) || fps < 0
    error('Tadpole:movie', 'fps must be a number > 0 (make a video) or 0 (images only).');
end

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
    if isempty(sediment)
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

% the video file: MPEG-4 if this platform has it, otherwise Motion JPEG AVI
if fps > 0
    videofile = [basename '.mp4'];
    try
        v = VideoWriter(videofile, 'MPEG-4');
    catch
        videofile = [basename '.avi'];
        warning('Tadpole:movie', 'MPEG-4 is not available on this platform; writing %s (Motion JPEG AVI) instead.', videofile);
        v = VideoWriter(videofile, 'Motion JPEG AVI');
    end
    v.Quality = 95; % 0-100: higher is sharper and larger
    v.FrameRate = fps;
    open(v);
    closeVideo = onCleanup(@() close(v)); % close the file even if a frame fails
end

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
    imagefile = [framefile num2str(framenum, '%04d') '.png'];
    print(p.fighandle, '-dpng', imagefile); % png

    % add the image to the video
    if fps > 0
        im = imread(imagefile);
        if size(im,3) == 1
            im = repmat(im, [1 1 3]); % grayscale to RGB
        end
        if framenum == 1
            framesize = size(im);
        elseif ~isequal(size(im), framesize)
            % the figure changed size (for example, the window was resized):
            % crop the image, or pad it with white, to the size of the first
            padded = 255*ones(framesize, 'like', im);
            nr = min(size(im,1), framesize(1));
            nc = min(size(im,2), framesize(2));
            padded(1:nr,1:nc,:) = im(1:nr,1:nc,:);
            im = padded;
        end
        % H.264 needs an even number of rows and columns: repeat the last
        % row or column of an odd-sized image
        if mod(size(im,1),2), im(end+1,:,:) = im(end,:,:); end
        if mod(size(im,2),2), im(:,end+1,:) = im(:,end,:); end
        writeVideo(v, im);
    end

end

% close the figure
close(p.fighandle);

% close the video file
if fps > 0
    clear closeVideo
    fprintf('Wrote %s: %d frames at %g fps (%.1f s).\n', videofile, framenum, fps, framenum/fps);
end
