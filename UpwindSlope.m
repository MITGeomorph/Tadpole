function S = UpwindSlope(p,g)

% calculate slope using a weighted average of the slopes to the neighbors
% specified in W

dr = [0 1 1 1 0 -1 -1 -1];
dc = [-1 -1 0 1 1 1 0 -1];
diag = sqrt(p.dx^2 + p.dy^2);
invd = 1./[p.dx diag p.dy diag p.dx diag p.dy diag];

% slopes to 8 nearest neighbors
S = zeros(p.K,p.J,8);
for k = 1:8
    S(:,:,k) = invd(k) * ( g.U - circshift(g.U,[dr(k) dc(k)]) );
end

% % assign zeros for fixed or mirror boundaries
% if bdy.left == 0 || bdy.left == 1
%     idx(:,1,1:3) = 0;
% end
% 
% if bdy.right == 0 || bdy.right == 1
%     idx(:,nc,7:9) = 0;
% end
% 
% if bdy.upper == 0 || bdy.upper == 1
%     idx(1,:,[1 4 7]) = 0;
% end
% 
% if bdy.lower == 0 || bdy.lower == 1
%     idx(nr,:,[3 6 9]) = 0;
% end

% multiply by weights and sum
S = sum(S.*g.W,3);


% % Optionally, use centered differencing on hillslopes
% Scentered = mexSlope(g.U,p.dx,p.dy,p.bvec);
% 
% S(~g.Channels) = Scentered(~g.Channels);
