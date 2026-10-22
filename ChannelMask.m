function Channels = ChannelMask(p, g)

% Channels = ChannelMask(p, g)
%
% Cells that count as channels: where rivers erode with p.channelThreshold =
% 'area', and the channel sinks of hillslope transport when channels do not
% diffuse (p.doChannelDiffusion = 0).
%
%   'area'   (the original rule)  A > p.Ac
%   'smooth' (default)            Omega > Omega_c, with Omega = (P*A)^m S^n and the
%            threshold stream power Omega_c = Ecr/Kr, or with sediment the cover-
%            weighted exp(-H/H*) Ecr/Kr + (1 - exp(-H/H*)) Ecs/Ks.
%
% g must hold the routing of DrainageArea (A, Q, W, Flooded, Zf). Slopes are
% taken on g.Z, with lake cells at the lake surface as in SPEqn; after a
% river step the weights are those it used (routed before the step's erosion).

if strcmp(p.channelThreshold, 'area')
    Channels = g.A > p.Ac;
    return
end

gs = g;
gs.Z(g.Flooded > 0) = g.Zf(g.Flooded > 0);
S = UpwindSlope(p, gs);
S(S < 0) = 0;
Om = g.Q.^p.m .* S.^p.n;
Omc = p.Ecr./p.Kr;
if p.doSediment
    e = exp(-g.H/p.Hstar);
    Omc = e.*Omc + (1 - e).*p.Ecs./p.Ks;
end
Channels = Om > Omc & g.C == 1;

end
