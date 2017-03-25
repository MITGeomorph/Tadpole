function [RHS Channels Courant] = NKWE(g,p)

% evaluates the RHS of the nonlinear kinematic wave term


% calculate drainage area
g = DrainageArea(p,g);

% calculate slope
S = UpwindSlope(p,g);
S(S<0) = 0; % This will only happen where we've routed flow uphill to flood depressions


% calculate Keff, which accounts for channel width < grid spacing
if p.wexp == 0
    Keff = p.Kf * p.Kw / p.dx;
elseif p.wexp == 1
    Keff = p.Kf * p.Kw * g.A / p.dx;
else
    Keff = p.Kf * p.Kw * g.A.^p.wexp ./ p.dx;
end


% calculate A^m*S^w
if p.m ~= 1 % otherwise we don't need to raise area to a power, saving time
    Am = g.A.^p.m;
else
    Am = g.A;
end

if p.w ~= 1 % ditto for slope   
    Sw = S.^p.w;
else
    Sw = S;
end

AmSw = Am.*Sw;


% record locations of channels (points that exceed threshold)
Channels = AmSw > p.thetac;



% calculate right-hand side of NKWE 

if p.thetac == 0
        
    RHS = -Keff .* AmSw;    
    
else
    
    % Case 1: if erosion goes as excess shear stress (conventional):
    RHS = -Keff .* (AmSw - p.thetac);
    
%     % Case 2: if erosion goes as shear stress:
%     RHS = -Keff .* AmSw;

    
    % In either case, this term is only nonzero where AmSw > thetac (bedrock can only erode)
    RHS = RHS.*(RHS < 0); 
    
end



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
% the time step, not the grid identified locally, since this call to NKWE
% could be one of the sub-steps in the RK2 method.

% if p.flood==1
%     RHS = RHS.*(~fl & g.C==1 & g.Channels); 
% else
%     RHS = RHS.*(~min & g.C==1 & g.Channels);
% end

if p.flood==1
    RHS = RHS.*(~g.fl & g.C==1 & Channels); 
else
    RHS = RHS.*(~g.minima & g.C==1 & Channels);
end



% maximum Courant #, neglecting threshold:
% if dz/dt = -Keff*A^m*S^n
% then C = Keff*A^m*S^(n-1)*dt/dx

Courant = max(Keff(:).*AmSw(:)./S(:))*sqrt(2)*p.dt/p.dx;
