function g = DrainageArea(p,g)

% compute upslope contributing area and upstream-averaged precipitation rate using the requested flow routing method
%
% mexD8, mexDinf and mexDms all take the same arguments. The optional 7th
% argument is the slope exponent for multiple-direction flow routing (used
% only by mexDms; flow is divided among downslope neighbors in proportion to
% slope^mfdexp). It is set by p.mfdexp if that field exists, and is 1 otherwise.
%
% g.Stack lists the linear indices of all cells, each before all of its
% donors; g.Receiver(i,j,k) is the linear index of the cell that receives
% the fraction g.W(i,j,k) of the flow from cell (i,j), or 0 if g.W(i,j,k) is 0.
% g.Zf is the surface after Priority-Flood filling (equal to g.Z where nothing
% was filled). g.Lake labels each flooded cell with the linear index of the
% spill cell of its depression, so all cells of one lake share a label; it is
% 0 for cells that were not flooded.
if isfield(p,'mfdexp')
    mfdexp = p.mfdexp;
else
    mfdexp = 1;
end
[g.A, g.Minima, g.Flooded, g.W, g.Pav, g.Stack, g.Receiver, g.Zf, g.Lake] = feval(['mex' p.routing], g.Z, p.dy/p.dx, 1*(~g.C), p.bvec, p.flood, p.P, mfdexp);

% multiply by cell area
g.A = p.dx*p.dy*g.A;

% water discharge Q = g.Pav .* g.A, assuming steady runoff with no losses
g.Q = g.Pav .* g.A;

% Test to make sure Stack is working: rebuild A and Q by looping through
% the inverted stack, and compare with the mex function. This loop is slow,
% so it is off by default.
checkStack = false;
if checkStack
    N = p.K*p.J;
    A = p.dx*p.dy*ones(p.K,p.J); % contribution from each cell to itself
    Q = p.dx*p.dy*p.P; % contribution from each cell to itself

    for t = N:-1:1
        c = g.Stack(t); % linear index of this cell
        for k = 1:8 % loop through neighbor directions
            e = c + N*(k-1); % linear index of direction k of this cell in W and Receiver
            if g.W(e) % if this direction has a nonzero weight
                r = g.Receiver(e); % linear index of this receiver neighbor
                A(r) = A(r) + g.W(e)*A(c); % Add the fraction of c's area that drains to r to A(r)
                Q(r) = Q(r) + g.W(e)*Q(c); % Same for water discharge
            end
        end
    end

    % On mirror boundaries, the mex functions double the flow that a boundary
    % cell receives from the interior (the mirror image of the donor also
    % drains to it), which this loop does not, so only compare without them.
    if ~any(p.bvec == 1)
        relErr = max(abs(A(:) - g.A(:)) ./ g.A(:));
        if relErr > 1e-9
            warning('DrainageArea:stack', 'Stack accumulation differs from mex drainage area (max relative error %g).', relErr);
        end
    end
end
