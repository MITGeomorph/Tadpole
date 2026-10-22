function g = Uplift(p,g)

% Uplift/subsidence relative to boundaries. Note that uplift occurs only where C == 1
g.Z = g.Z + p.dt.*p.U.*g.C;
