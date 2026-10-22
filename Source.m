function g = Source(p,g)

% Source terms or perturbations (Uplift/subsidence, etc.)

% surface uplift relative to boundaries
g = Uplift(p,g);

% uplift raises the bedrock surface by the same amount; sediment thickness is unchanged
if p.doSediment
    g.R = g.R + p.dt.*p.U.*g.C;
end
