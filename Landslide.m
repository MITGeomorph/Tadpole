function g = Landslide(p,g)

sumdz = 1;

while sumdz

    [sorted order] = sort(g.U(:));
    [g.U dz] = mexLandslide(g.U,order,p.dx,p.dy,p.Sc,p.bvec,1*(~g.C));
    sumdz = nansum(dz(:));
    
end

% we could possibly avoid running more than once if there are no local
% minima. Check to see how long it's taking, and see if it makes a
% difference. In practice, it seems that there can be aftereffects even
% once all local minima have failed