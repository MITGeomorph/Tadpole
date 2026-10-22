function DrawPlotMovie(iter,p,g)

switch p.plottype

    case 'elevation' % map view with elevation as color
        
        figure(p.fighandle)
        imagesc(g.x,g.y,g.Z,[p.zmin,p.zmax])
        title(['Iteration=' num2str(iter) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        colorbar

        set(gca,'visible','off')
        colorbar off
        
        drawnow

    case {'hillshade','sediment'} % Top-down shaded relief, colored by elevation ('hillshade') or by
                                  % sediment thickness ('sediment')

        figure(p.fighandle)

        v = 5; % vertical exaggeration (scalar)
        az = 315; % light azimuth in degrees CCW from +x

        if strcmp(p.plottype,'sediment')
            h = surf(g.x, g.y, v*g.Z, g.H);   % color by sediment thickness
            lims = [p.Hmin p.Hmax];
        else
            h = surf(g.x, g.y, v*g.Z, g.Z);   % color by original elevation
            lims = [p.zmin p.zmax];
        end
        set(h, 'EdgeColor', 'none');
        
        title(['Iteration=' num2str(iter) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])

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

        figure(p.fighandle)
        PerspectiveView(p, g, p.persp, ['t = ' num2str(round(p.t)) ' yr']);

end
