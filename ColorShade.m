function [rgb,hi] = ColorShade(Z,C,az,alt,zfac,doplot,cmin,cmax)

% ColorShade(Z,C,az,alt,zfac,doplot)
%
% Calculates, and optionally plots, a colored shaded relief map using the
% current colormap
%
%   Arguments:
%        Z.grid     used to calculate shaded relief
%        C.grid     used to calculate colors (defaults to Z)
%        az         direction of light source in deg E of N (default 315)
%        alt        altitude of the lighting source in degrees above 
%                   horizontal (default 45)
%        zfac       elevation scaling z-factor (default 1)
%        doplot     specifies whether to plot the hillshade (default 1)
%
%        hi         handle of the image object containing the map
%        rgb        the rgb values for the colored shaded relief map
%        cmin,cmax  min and max values for the color scale, in units of C

% Based on imshade.m by Kelsey Jordahl, Marymount Manhattan College, 2009
% Copyright (c) 2009, Kelsey Jordahl
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Modifications Copyright (C) 2010 Taylor Perron <perron@mit.edu>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation. You should have received a copy of the GNU 
% General Public License along with this program.  If not, see 
% http://www.gnu.org/licenses.

% Set default parameters if not specified
if nargin < 6
    doplot = 1;
    if nargin < 5
        zfac = 1;
        if nargin < 4 
            alt = 45;
            if nargin < 3
                az = 315;
                if nargin < 2
                    C = Z;
                end
            end
        end
    end
end

if nargin < 8
    
    cmin = min(C.grid(:));
    cmax = max(C.grid(:));
    
end

nc = length(colormap);

% get rgb values for color grid

% cmin = min(C.grid(:));
% cmax = max(C.grid(:));
% C.grid=C.grid-min(C.grid(:));
% C.grid=C.grid*255/max(C.grid(:));
% rgb=ind2rgb(uint8(C.grid*length(colormap)/255),colormap);
C.grid=(C.grid-cmin)*255/(cmax-cmin);
C.grid(C.grid<0)=0;
C.grid(C.grid>255)=255;
rgb=ind2rgb(uint8(C.grid*nc/255),colormap);


% calculate hillshade and scale to values between 0 and 1
shade = Shade(Z,az,alt,zfac,0); 
shade = shade/255;

% convert to hsv and set value to scaled hillshade, then convert back to
% rgb
hsv=rgb2hsv(rgb);
hsv(:,:,3) = shade;
rgb = hsv2rgb(hsv);

% draw a plot if requested, and add a colorbar with appropriate labels
if doplot
    hi = image(Z.x,Z.y,rgb);
    axis image; axis xy
    h = colorbar;
    t = get(h,'ytick');
    t = cmin + t/nc*(cmax-cmin);
    set(h,'yticklabel',t);
end
