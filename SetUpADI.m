function g = SetUpADI(p,g)

N=p.J*p.K;

% Create vectors of row, column, and values
if p.doChannelDiffusion
    dodiff = 1*g.C;
else
    dodiff = 1*(g.C & ~g.Channels);
end
[Ar Ac Av Br Bc Bv amax bmax] = mexGetMats(p.D,p.dt,p.dx,p.dy,p.bvec,dodiff);

% eliminate unused elements
Ar = Ar(1:amax);
Ac = Ac(1:amax);
Av = Av(1:amax);

Br = Br(1:bmax);
Bc = Bc(1:bmax);
Bv = Bv(1:bmax);


% Create the sparse matrices used in the Crank-Nicolson ADI steps
A = sparse(Ar,Ac,Av,N,N);
B = sparse(Br,Bc,Bv,N,N);

I = sparse(1:N,1:N,1);

g.Aplus = I + A;
g.Aminus = I - A;
g.Bplus = I + B;
g.Bminus = I - B;
