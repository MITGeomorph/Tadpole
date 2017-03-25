function DrawPlot(n,p,g)

switch p.plottype

    case 'drainage area' % map view with log(A) as color
        
        figure(p.fighandle)
        imagesc(g.x,g.y,log10(eval(['mex' p.routing '(g.U,p.dy/p.dx,1*(~g.C),p.bvec,p.flood)'])))
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow                

    case 'curvature' % map view with normalized curvature as color
        
        figure(p.fighandle)
        imagesc(g.x,g.y,curvn(g.U))
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow                

    case 'elevation' % map view with elevation as color
        
        figure(p.fighandle)
        imagesc(g.x,g.y,g.U)
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        colorbar

%         set(gca,'visible','off')
%         colorbar off
        
        drawnow
        
    case 'contour' % normalized contour map
        
        figure(p.fighandle)
        ncont=10;
        contour(g.X,g.Y,(g.U-min(g.U(:)))/range(g.U(:)),ncont)
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow        

    case 'shade' % shaded relief
        
        figure(p.fighandle)
        Z.grid = flipud(g.U);
        Z.x = g.x;
        Z.y = g.x;
        Z.dx = p.dx;
        Z.dy = p.dy;
        az = 315;
        alt = 45;
        Shade(Z,az,alt,p.zfactor);
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow

    case 'color shade' % colored shaded relief
        
        figure(p.fighandle)
        Z.grid = flipud(g.U);
        Z.x = g.x;
        Z.y = g.x;
        Z.dx = p.dx;
        Z.dy = p.dy;
        az = 315;
        alt = 45;
%         ColorShade(Z,Z,az,alt,p.zfactor);
        [~,h] = ColorShade(Z,Z,az,alt,p.zfactor);
        A = zeros(size(p.F));
        A(p.F==1) = 0;
        A(p.F==0) = 1;
        set(h,'alphadata',flipud(A))
%         colorbar off
%         set(gca,'visible','off')
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow

    case 'mesh' % mesh plot

        figure(p.fighandle)
        surf(g.X,g.Y,flipud(g.U)); 
        vertexag = [1 max(g.Y(:))/max(g.X(:)) p.zfactor*range(g.U(:))/max(g.X(:))]; 
        set(gca, 'PlotBoxAspectRatio', vertexag);
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow
 

    otherwise % case 1 and default: perspective surface plot

        figure(p.fighandle)
        surf(g.X,g.Y,flipud(g.U)); 
        shading interp
        camlight left 
        lighting phong
        vertexag = [1 max(g.Y(:))/max(g.X(:)) p.zfactor*range(g.U(:))/max(g.X(:))]; 
        set(gca, 'PlotBoxAspectRatio', vertexag);
        title(['Iteration=' num2str(n) ', t=' num2str(round(p.t)) ', dt=' num2str(round(p.dt))])
        drawnow

end