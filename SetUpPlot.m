function [p g] = SetUpPlot(p,g)

p.fighandle = figure; 
figure(p.fighandle)

switch p.plottype

    case 'drainage area' % map view with log(A) as color
        
        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        imagesc(g.x,g.y,log10(p.dx*p.dy*eval(['mex' p.routing '(g.U,p.dy/p.dx,1*(~g.C),p.bvec,p.flood)'])))
        axis image
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        drawnow                

    case 'curvature' % map view with normalized curvature as color
        
        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        imagesc(g.x,g.y,curvn(g.U))
        axis image
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        drawnow
        
    case 'elevation' % map view with elevation as color
        
        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        imagesc(g.x,g.y,g.U)
        axis image
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        colorbar
        
%         set(gca,'visible','off')
%         colorbar off

        drawnow        

    case 'contour' % normalized contour map
        
        ncont=10;
        [g.X g.Y] = meshgrid(p.dx*(0:p.J-1),p.dy*(0:p.K-1));        
        contour(g.X,g.Y,(g.U-min(g.U(:)))/range(g.U(:)),ncont)
        axis image
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        drawnow        
        
    case 'shade' % shaded relief

        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        Z.grid = flipud(g.U);
        Z.x = g.x;
        Z.y = g.y;
        Z.dx = p.dx;
        Z.dy = p.dy;
        az = 315;
        alt = 45;
        if ~isfield(p,'zfactor')
            p.zfactor = 1;
        end
        Shade(Z,az,alt,p.zfactor);
        axis image

%         hold on
%         ind = find(p.F);
%         [i j] = ind2sub([p.K p.J],ind);
%         plot(j,i,'.w')
%         set(gca,'ydir','reverse')
%         hold off
        
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1]);
        drawnow        

    case 'color shade' % colored shaded relief

        g.x = p.dx*(0:p.J-1);
        g.y = p.dy*(0:p.K-1);
        Z.grid = flipud(g.U);
        Z.x = g.x;
        Z.y = g.y;
        Z.dx = p.dx;
        Z.dy = p.dy;
        az = 315;
        alt = 45;
        if ~isfield(p,'zfactor')
            p.zfactor = 1;
        end
%         ColorShade(Z,Z,az,alt,p.zfactor);
        [~,h] = ColorShade(Z,Z,az,alt,p.zfactor);
        A = zeros(size(p.F));
        A(p.F==1) = 0;
        A(p.F==0) = 1;
        set(h,'alphadata',flipud(A))

%         set(gca,'visible','off')
%         colorbar off

        axis image

%         hold on
%         ind = find(p.F);
%         [i j] = ind2sub([p.K p.J],ind);
%         plot(j,i,'.w')
%         set(gca,'ydir','reverse')
%         hold off
        
        title(['Iteration= 0, t=' num2str(p.t)])
        set(gca, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        drawnow        

    case 'mesh' % mesh plot

        [g.X g.Y] = meshgrid(p.dx*(0:p.J-1),p.dy*(0:p.K-1));
        surf(g.X,g.Y,flipud(g.U));
        title(['Iteration= 0, t=' num2str(p.t)])
        if ~isfield(p,'zfactor')
            p.zfactor = 1;
        end
        vertexag = [1 max(g.Y(:))/max(g.X(:)) p.zfactor*range(g.U(:))/max(g.X(:))]; 
        set(gca, 'PlotBoxAspectRatio', vertexag, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        drawnow

    otherwise % perspective surface plot

        [g.X g.Y] = meshgrid(p.dx*(0:p.J-1),p.dy*(0:p.K-1));
        surf(g.X,g.Y,flipud(g.U));
        title(['Iteration= 0, t=' num2str(p.t)])
        shading interp
        camlight left
        lighting phong
        if ~isfield(p,'zfactor')
            p.zfactor = 1;
        end
        vertexag = [1 max(g.Y(:))/max(g.X(:)) p.zfactor*range(g.U(:))/max(g.X(:))]; 
        set(gca, 'PlotBoxAspectRatio', vertexag, 'nextplot', 'replacechildren');
        set(p.fighandle, 'color', [1 1 1], 'colormap', jet);
        drawnow
    
end
