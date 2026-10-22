function [RHS, g] = SPEqn(g,p)

% evaluates the RHS of the stream power equation


% calculate drainage area
g = DrainageArea(p,g);

% calculate slope. A river flowing into a lake erodes down to the lake's
% water surface (the filled surface), not to the lake floor, because stream
% power depends on the water surface slope. (With p.flood = 0 there are no
% lakes and g.Zf = g.Z.)
gs = g;
gs.Z(g.Flooded > 0) = g.Zf(g.Flooded > 0);
S = UpwindSlope(p,gs);
S(S<0) = 0; % This will only happen where we've routed flow uphill to flood depressions


Keff = p.Kr; % erodibility (bedrock; the explicit solver does not handle sediment)

% calculate (P*A)^m*S^n
if p.m ~= 1 % otherwise we don't need to raise area to a power, saving time
    PAm = (g.Pav.*g.A).^p.m;
else
    PAm = g.Pav.*g.A;
end

if p.n ~= 1 % ditto for slope   
    Sn = S.^p.n;
else
    Sn = S;
end

PAmSn = PAm.*Sn;

% record locations of channels: A > p.Ac with p.channelThreshold = 'area'; with
% 'smooth' (default), Omega > Ecr/Kr (ChannelMask)
g.Channels = ChannelMask(p, g);
if strcmp(p.channelThreshold, 'area')
    active = g.Channels;
else
    active = true(size(g.Z)); % every cell follows the erosion law with the threshold
end
Ecr = ErosionThresholds(p);

% calculate right-hand side of the stream power equation: erosion rate with the smoothed threshold
% Ecr, F(Kr Omega; Ecr) (= Kr Omega without a threshold)

RHS = -SmoothThreshold(Keff .* PAmSn, Ecr);    

% RHS = RHS.*(RHS < 0); % This term cannot be positive (bedrock can only erode)
    
% no-erosion rules 

% no erosion in:
%
% 1. areas that are flooded (or local minima, if the user did not request
% to flood depressions when doing flow routing). We assume these are areas of
% temporary deposition.
%
% 2. Fixed boundaries and other locations requested by user, as recorded in the matrix C.
%
% 3. Non-channels. Note that if we use a midpoint method like RK2, we 
% should use the channel grid determined at the beginning of
% the time step, not the grid identified locally, since this call to SPEqn
% could be one of the sub-steps in the RK2 method.
% 
% 4. Also, RHS cannot be positive (bedrock can only erode)

if p.flood==1
    RHS = RHS.*(~g.Flooded & g.C==1 & active & RHS < 0); 
else
    RHS = RHS.*(~g.Minima & g.C==1 & active & RHS < 0);
end



% maximum Courant #, neglecting threshold:
% if dz/dt = -Keff*A^m*S^n
% then C = Keff*A^m*S^(n-1)*dt/dx

% Cells with S = 0 do not erode and have no erosion wave: they contribute 0
% (0/0 gave NaN there, which max ignored)
Kv = Keff.*ones(size(S));
pos = S > 0;
g.Courant = max([0; Kv(pos).*PAmSn(pos)./S(pos)])*sqrt(2)*p.dt/p.dx;
