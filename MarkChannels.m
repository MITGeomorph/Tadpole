function [p g] = MarkChannels(p,g)

% [p g] = MarkChannels(p,g)
%
% Identifies the locations of channels based on fluvial incision threshold,
% and records them in a matrix, g.Channels

if p.thetac > 0
    
    % calculate drainage area
    g = DrainageArea(p,g);
    
    % calculate slope
    S = UpwindSlope(p,g);
    S(S<0) = 0; % This will only happen where we've routed flow uphill to flood depressions
    
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
    g.Channels = AmSw > p.thetac;
    
else
    
    % If there is no fluvial incision threshold, channel incision can occur everywhere
    g.Channels = ones(p.K,p.J);
    
end
