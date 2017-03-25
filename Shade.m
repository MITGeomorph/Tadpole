function h = Shade(Z,az,alt,zfac,doplot)

% h = Shade(Z,az,alt,zfac,doplot)
%
% Calculate hillshade for a digital elevation model (DEM)
%
% where: Z.grid     DEM to calculate hillshade for
%        Z.dx, Z.dy grid spacings
%        az         direction of light source in deg E of N (default 315)
%        alt        altitude of the lighting source in degrees above 
%                   horizontal (default 45)
%        zfac       elevation scaling z-factor (default 1)
%        doplot     specifies whether to plot the hillshade (default 1)
%

% Based on hillshade.m by Felix Hebeler, Dept. of Geography, University 
% of Zurich, February 2007.
% Modified by Andrew Stevens (astevens@usgs.gov), 5/04/2007
% Modified by Taylor Perron (perron@mit.edu), July 2010

% Set default parameters if not specified
if nargin < 5
    doplot = 1;
    if nargin < 4
        zfac = 1;
        if nargin < 3 
            alt = 45;
            if nargin < 2
                az = 315;
            end
        end
    end
end

% convert light azimuth to polar angle in radians
az = 360.0 - az + 90;  
az(az>=360) = az - 360;
az = az * (pi/180);

% convert lighting altitude to zenith angle in radians
alt = (90-alt) * (pi/180);

% calculate slope and aspect in radians
[fx,fy] = gradient(Z.grid,Z.dx,Z.dy); 
[asp,grad] = cart2pol(fy,fx);
grad = atan(zfac*grad); % scale slope by z-factor

% convert aspect angle
asp(asp<pi) = asp(asp<pi)+(pi/2);
asp(asp<0) = asp(asp<0)+(2*pi);

% hillshade calculation
h = 255.0*( (cos(alt).*cos(grad) ) + ( sin(alt).*sin(grad).*cos(az-asp)) ); % ESRIs algorithm
h(h<0) = 0; % set hillshade values to min of 0.

h = setborder(h,1,NaN); % set border cells to NaN

if doplot

    imagesc(Z.x,Z.y,h);
    colormap gray
    axis image; axis xy
    
end



% -- Subfunction --
function grid = setborder(grid,bs,bvalue)
grid(1:bs,:)=bvalue; %toprows
grid(size(grid,1)-bs+1:size(grid,1),:)=bvalue; %bottom rows
grid(:,1:bs)=bvalue; %left cols
grid(:,size(grid,2)-bs+1:size(grid,2))=bvalue;