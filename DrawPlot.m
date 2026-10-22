function DrawPlot(iter,p,g)

switch p.plottype

    case 'elevation' % map view with elevation as color
        
        figure(p.fighandle)
        imagesc(g.x,g.y,g.Z)
        title(['Iteration=' num2str(iter) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        colorbar
        
        drawnow

    case {'hillshade','sediment'} % Top-down shaded relief, colored by elevation ('hillshade') or by
                                  % sediment thickness ('sediment')

        figure(p.fighandle)

        v = 5; % vertical exaggeration (scalar)
        az = 315; % light azimuth in degrees CCW from +x

        if strcmp(p.plottype,'sediment')
            h = surf(g.x, g.y, v*g.Z, g.H);   % color by sediment thickness
        else
            h = surf(g.x, g.y, v*g.Z, g.Z);   % color by original elevation
        end
        set(h, 'EdgeColor', 'none');
        
        title(['Iteration=' num2str(iter) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])

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

        figure(p.fighandle)
        PerspectiveView(p, g, p.persp, ['Iteration=' num2str(iter) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))]);

end