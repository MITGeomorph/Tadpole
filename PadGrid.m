function M = PadGrid(M, p)

% Pads a K x J grid with one row or column on each side according to the
% boundary conditions in p.bvec = [left right upper lower]:
% 0 = fixed (padded with zeros), 1 = mirror (reflected about the edge cell),
% 2 = periodic (wrapped from the opposite edge).

[K, J] = size(M);
Lpad = zeros(K,1);
Rpad = zeros(K,1);
Tpad = zeros(1,J+2); % The +2 is because we pad L and R boundaries first.
Bpad = zeros(1,J+2);

if p.bvec(1) == 1, Lpad = M(:,2); elseif p.bvec(1) == 2, Lpad = M(:,J); end
if p.bvec(2) == 1, Rpad = M(:,J-1); elseif p.bvec(2) == 2, Rpad = M(:,1); end
M = [Lpad M Rpad];

if p.bvec(3) == 1, Tpad = M(2,:); elseif p.bvec(3) == 2, Tpad = M(K,:); end
if p.bvec(4) == 1, Bpad = M(K-1,:); elseif p.bvec(4) == 2, Bpad = M(1,:); end
M = [Tpad; M; Bpad];
