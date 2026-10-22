function S = UpwindSlope(p,g)

% calculate slope using a weighted average of the slopes to the neighbors
% specified in W
%
% Uses g.Z, g.W, g.Flooded, p.dx, p.dy, p.K, p.J and p.bvec. S is set to zero
% in flooded cells (where g.Flooded is nonzero, as returned by the flow routing
% functions in DrainageArea.m), because they do not erode.
%
% On mirror boundaries, the flow routing functions halve the weights of
% flow from a boundary cell into the interior (FixMirrorWeights), so that
% drainage area is counted correctly. Those weights are doubled again here,
% so the slope of a mirror boundary cell is a weighted average of its
% slopes with weights that sum to 1, as everywhere else. Corner cells are
% included; at a corner between two mirror boundaries, a diagonal weight into
% the interior was halved twice, so it is doubled twice here.

dr = [0 1 1 1 0 -1 -1 -1];
dc = [-1 -1 0 1 1 1 0 -1];
diag = sqrt(p.dx^2 + p.dy^2);
invd = 1./[p.dx diag p.dy diag p.dx diag p.dy diag];

% directions (1-8) whose weights are halved on each mirror boundary, as in
% FixMirrorWeights: [left; right; upper; lower]
halved = [8 1 2; 4 5 6; 6 7 8; 2 3 4];
mirror = p.bvec(:)' == 1;

% sum the weighted slopes to the 8 nearest neighbors, one direction at a time
S = zeros(p.K,p.J);
for k = 1:8
    Wk = g.W(:,:,k);
    if mirror(1) && any(halved(1,:) == k), Wk(:,1)   = 2*Wk(:,1);   end
    if mirror(2) && any(halved(2,:) == k), Wk(:,end) = 2*Wk(:,end); end
    if mirror(3) && any(halved(3,:) == k), Wk(1,:)   = 2*Wk(1,:);   end
    if mirror(4) && any(halved(4,:) == k), Wk(end,:) = 2*Wk(end,:); end
    % circshift wraps around every edge, but weights across fixed and mirror
    % edges are zero, so only the periodic wrap-around contributes
    S = S + Wk .* (invd(k) * (g.Z - circshift(g.Z,[dr(k) dc(k)])));
end

% flooded cells do not erode, so their slope is set to zero
if isfield(g,'Flooded') && ~isempty(g.Flooded)
    S(g.Flooded ~= 0) = 0;
end


% % Optionally, use centered differencing on hillslopes
% Scentered = mexSlope(g.Z,p.dx,p.dy,p.bvec);
%
% S(~g.Channels) = Scentered(~g.Channels);
