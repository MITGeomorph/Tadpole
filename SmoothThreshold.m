function F = SmoothThreshold(x, Ec)

% F = SmoothThreshold(x, Ec)
%
% Erosion rate with a smoothed threshold (Shobe et al., 2017; Landlab's
% SPACE):
%
%   F(x; Ec) = x - Ec*(1 - exp(-x/Ec))
%
% where x = K*Omega is the erosion rate without a threshold and Ec (m/yr) the
% threshold erosion rate. F ~ x^2/(2 Ec) for x << Ec, F ~ x - Ec for x >> Ec,
% and F = x for Ec = 0. Both F and its derivative are continuous, so there is
% no sharp switch between eroding and non-eroding cells. For x/Ec < 1e-3 a
% series is used (the direct formula loses accuracy there); mexSPACE.c uses
% the same formula, so the two can be compared to round-off.

if Ec <= 0
    F = x;
    return
end
y = x/Ec;
F = x - Ec*(1 - exp(-y));
small = y < 1e-3;
ys = y(small);
F(small) = Ec*ys.^2.*(0.5 - ys/6 + ys.^2/24);

end
