function h = PerspectiveView(p, g, h, ttl)

% h = PerspectiveView(p, g, h, ttl)
%
% Oblique 3-D view of the topography for p.plottype = 'perspective': smooth
% shading, two lights, an earth-tone colormap and a dark background. With h
% empty it builds the view in the current axes; otherwise it updates the
% surface in h (much faster than redrawing).
%
% Vertical exaggeration: p.plotVE if set; otherwise chosen each time so that
% the relief spans 1/8 of the larger horizontal dimension (at least 1), and
% shown in the title. Color limits: [p.zmin p.zmax] if set (TadpoleMovie
% fixes them for all frames), otherwise the current range of elevations.
% Grids larger than ~600 cells on a side are thinned for display only.

Z = g.Z;
[K, J] = size(Z);
s = max(1, ceil(max(K, J)/600)); % display stride
ii = 1:s:K; jj = 1:s:J;
x = p.dx*(jj - 1); y = p.dy*(ii - 1);
Zs = Z(ii, jj);

% vertical exaggeration
relief = max(Zs(:)) - min(Zs(:));
if isfield(p,'plotVE') && ~isempty(p.plotVE)
    ve = p.plotVE;
else
    ve = max(1, max(p.dx*(J-1), p.dy*(K-1))/8/max(relief, eps));
end

if isfield(p,'zmin') && isfield(p,'zmax') && p.zmax > p.zmin
    lims = [p.zmin p.zmax];
else
    lims = [min(Zs(:)) max(Zs(:)) + eps];
end

if isempty(h)
    ax = gca; cla(ax);
    bg = [0.32 0.34 0.43]; % ParaView's default background
    set(gcf, 'Color', bg);
    h.land = surf(ax, x, y, ve*Zs, Zs, 'EdgeColor', 'none', 'FaceColor', 'interp', ...
        'FaceLighting', 'gouraud', 'AmbientStrength', 0.35, 'DiffuseStrength', 0.8, ...
        'SpecularStrength', 0.1, 'SpecularExponent', 10, 'BackFaceLighting', 'lit');
    colormap(ax, EarthColors(256));
    axis(ax, 'equal', 'tight', 'off');
    axis(ax, 'vis3d');              % keep the aspect ratio when the view is rotated
    set([ax h.land], 'Clipping', 'off'); % zooming narrows the axis limits; draw the whole surface anyway
    set(ax, 'YDir', 'reverse');     % row 1 at the top, as in the map views
    view(ax, -35, 35);
    camproj(ax, 'perspective');     % perspective rather than orthographic projection
    camlight(ax, -60, 40);          % key light
    camlight(ax, 120, 25);          % fill light
    h.ax = ax;
    cb = colorbar(ax); cb.Color = [0.9 0.9 0.9]; cb.Label.String = 'elevation (m)';
else
    set(h.land, 'XData', x, 'YData', y, 'ZData', ve*Zs, 'CData', Zs);
end
clim(h.ax, lims);
zr = ve*[min(Zs(:)) max(Zs(:)) + eps];
zlim(h.ax, zr + [-0.02 0.02]*diff(zr));
title(h.ax, sprintf('%s   (vertical exaggeration %.3g)', ttl, ve), 'Color', [0.9 0.9 0.9]);
drawnow

end

function c = EarthColors(n)
% green lowlands, tan and brown slopes, grey-white summits
k = [0    0.20 0.42 0.25
     0.30 0.45 0.58 0.30
     0.55 0.78 0.72 0.50
     0.78 0.58 0.45 0.33
     1.00 0.95 0.95 0.95];
t = linspace(0, 1, n)';
c = [interp1(k(:,1), k(:,2), t) interp1(k(:,1), k(:,3), t) interp1(k(:,1), k(:,4), t)];
end
