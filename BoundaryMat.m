function [p g] = BoundaryMat(p,g)

% create a matrix C that is 1 where uplift and erosion should occur, 
% 0 where it shouldn't
%
% and a vector bvec representing boundary conditions. bvec is 
% [left right upper lower] 
% with
% 0 = fixed elevation
% 1 = mirror
% 2 = periodic


g.C=ones(p.K,p.J);
p.bvec=zeros(1,4);

switch p.bdy.left

    case 'fixed'
        g.C(:,1)=0;
        p.bvec(1)=0;
    case 'mirror'
        p.bvec(1)=1;
    case 'periodic'
        p.bvec(1)=2;
end

switch p.bdy.right

    case 'fixed'
        g.C(:,p.J)=0;
        p.bvec(2)=0;
    case 'mirror'
        p.bvec(2)=1;
    case 'periodic'
        p.bvec(2)=2;
end

switch p.bdy.upper

    case 'fixed'
        g.C(1,:)=0;
        p.bvec(3)=0;
    case 'mirror'
        p.bvec(3)=1;
    case 'periodic'
        p.bvec(3)=2;
end

switch p.bdy.lower

    case 'fixed'
        g.C(p.K,:)=0;
        p.bvec(4)=0;
    case 'mirror'
        p.bvec(4)=1;
    case 'periodic'
        p.bvec(4)=2;
end


% Set additional fixed points

if isfield(p,'F') % if user specified additional fixed points
    g.C(p.F==1) = 0;
end

