function g = ADI(p,g)

% Crank-Nicolson (2nd-order accurate and unconditionally stable), 
% using an alternating-direction implicit (ADI) scheme
    
g.U = g.Aminus\(g.Bplus*g.U(:));
g.U = g.Bminus\(g.Aplus*g.U(:));
g.U = reshape(g.U,p.K,p.J);
