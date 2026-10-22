function [p, g] = SetUpPlotMovie(p,g)

p.fighandle = figure; 
figure(p.fighandle)

switch p.plottype
        
    case 'elevation' % map view with elevation as color
        
        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);        
        imagesc(g.x,g.y,g.Z,[p.zmin,p.zmax])
        axis image
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', parula);
        colorbar
        
        set(gca,'visible','off')
        colorbar off

        drawnow        

    case {'hillshade','sediment'} % Top-down shaded relief, colored by elevation ('hillshade') or by
                                  % sediment thickness ('sediment')

        v = 5; % vertical exaggeration (scalar)
        az = 315; % light azimuth in degrees CCW from +x

        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        if strcmp(p.plottype,'sediment')
            h = surf(g.x, g.y, v*g.Z, g.H);   % color by sediment thickness
            lims = [p.Hmin p.Hmax];
        else
            h = surf(g.x, g.y, v*g.Z, g.Z);   % color by original elevation
            lims = [p.zmin p.zmax];
        end
        set(h, 'EdgeColor', 'none');

        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1]);
        
        colormap(parula)                % choose colormap as desired
        shading interp
        lighting gouraud
        material([0.6 0.8 0.0])         % [ambient diffuse specular]
        axis equal tight
        view(2)                         % top-down

        clim(lims);
        set(gca, 'CLim', lims);
        
        % place a distant light at az (CCW from +x) and 45° elevation
        camlight(az, 45)
        colorbar

        set(gca,'visible','off')
        colorbar off

        drawnow

    case 'perspective' % oblique 3-D view with lighting, colored by elevation (PerspectiveView)

        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        p.persp = PerspectiveView(p, g, [], ['t = ' num2str(round(p.t)) ' yr']);

end
