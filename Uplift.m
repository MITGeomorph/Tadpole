function g = Uplift(p,g)

% Uplift/subsidence relative to boundaries. Note that uplift occurs only where C == 1
g.U = g.U + p.dt.*p.E.*g.C;
