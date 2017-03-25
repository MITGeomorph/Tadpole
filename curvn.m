function C = curvn(M)

% C = curvn(M)
%  creates a map of normalized (relative to max and min) curvature,
%  periodic in the x direction

% % pad M in x direction
% M = [M(:,end) M M(:,1)];

% calculate discrete laplacian
C = 4*del2(M);

% % extract portion corresponding to original dimensions of M
% C = C(:,2:end-1);

% normalize curvature relative to max & min values
C = -(C<0).*(C/min(C(:))) + (C>=0).*(C/max(C(:)));

