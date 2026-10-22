function [p, g] = SetUpPlot(p,g)

% title: iteration 0 at the start of a run; at the end (Snapshot), the number of
% iterations done
n = 0;
if isfield(p,'iterations'), n = p.iterations; end
ttl = ['Iteration=' num2str(n) ', t=' num2str(round(p.t))];

p.fighandle = figure; 
figure(p.fighandle)

switch p.plottype
        
    case 'elevation' % map view with elevation as color
        
        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        imagesc(g.x,g.y,g.Z)
        axis image
        title(ttl)
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', parula);
        colorbar
        
        drawnow        

    case {'hillshade','sediment'} % Top-down shaded relief, colored by elevation ('hillshade') or by
                                  % sediment thickness ('sediment')

        v = 5; % vertical exaggeration (scalar)
        az = 315; % light azimuth in degrees CCW from +x

        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        if strcmp(p.plottype,'sediment')
            h = surf(g.x, g.y, v*g.Z, g.H);   % color by sediment thickness
        else
            h = surf(g.x, g.y, v*g.Z, g.Z);   % color by original elevation
        end
        set(h, 'EdgeColor', 'none');

        title(ttl)
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1]);
        
        colormap(parula)                % choose colormap as desired
        shading interp
        lighting gouraud
        material([0.6 0.8 0.0])         % [ambient diffuse specular]
        axis equal tight
        view(2)                         % top-down
        
        % place a distant light at az (CCW from +x) and 45° elevation
        camlight(az, 45)
        cb = colorbar;
        if strcmp(p.plottype,'sediment'), cb.Label.String = 'sediment thickness (m)'; end

        drawnow

    case 'perspective' % oblique 3-D view with lighting, colored by elevation (PerspectiveView)

        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        p.persp = PerspectiveView(p, g, [], ttl);

end
