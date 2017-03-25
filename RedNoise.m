function M = RedNoise(ny,nx,beta,variance,periodic)

% M = RedNoise(ny,nx,beta,variance)
%
%  Generates a self-affine surface with dimensions (ny,nx) characterized by
%  the relationship
%
%  P = f^-beta
%
%  where P is spectral power and f is spatial frequency. The fractal 
%  dimension D is related to beta by
%
%  D = (7 - beta)/2
%
%  topography typically falls in the range 1.5 > beta > 2.5. The surface
%  heights are scaled such that the total variance is variance. beta = 0 
%  creates white noise; more positive values concentrate more variance at 
%  longer wavelengths. If periodic is 1, the output will be periodic at the
%  edges (works best if ny == nx).


% argument checking
if nargin < 1
    error('At least one input argument is required')
elseif nargin < 2
    nx = ny;
    beta = 2;
    variance = 1;
elseif nargin < 3
    beta = 2;
    variance = 1;
elseif nargin < 4
    variance = 1;
elseif nargin < 5
    periodic = 0;
end

if periodic
    
    if rem(ny,2) || rem(nx,2)
        error('ny and nx must be even')
    end
    
    % grid size
    [x,y] = meshgrid(1:nx,1:ny);
    
    % Generate a random matrix of Fourier components. Note that if you need a
    % different starting surface each time you launch a Matlab session, you
    % should initialize the random number generator with something like:
    % s = sum(100*clock);
    % rand('seed',s);
    F = (rand(ny,nx)-0.5) + 1i* (rand(ny,nx)-0.5);
    
    % Identify the DC component
    xc = nx/2+1; yc = ny/2+1; % matrix indices of zero frequency
    
    % make a matrix with entries proportional to the frequency
    freq = sqrt( (x-xc).^2 + (y-yc).^2 );
    
    % Reduce high freq components by a power law f^-beta
    %to make a fractal surface
    F = F .* freq .^ -beta;
    
    % Set the DC level (= mean of the elevations) to zero
    F(yc,xc) = 0;
    
    % Take the inverse FFT
    M = real(ifft2(ifftshift(F)));
        
    % scale elevations to a specified total variance
    M = M*sqrt(variance)/std(M(:));
    
    % % alternatively, scale elevations to produce a given relief, such that range(M(:)) = relief
    % M = M*relief/range(M(:));
    
    % % Optional: detrend the surface by subtracting a least-squares plane
    % M = detrend(M);
    
else % user did not request periodic output
    
    % grid size
    n = 2.^(ceil(log(max([nx ny]))/log(2))); % next largest power of 2       
    [x,y] = meshgrid(1:n);

    % Generate a random matrix of Fourier components. Note that if you need a 
    % different starting surface each time you launch a Matlab session, you 
    % should initialize the random number generator with something like:
    % s = sum(100*clock);
    % rand('seed',s);
    F = (rand(n)-0.5) + 1i* (rand(n)-0.5);  

    % Identify the DC component
    nc = n/2+1;

    % make a matrix with entries proportional to the frequency
    freq = sqrt( (x-nc).^2 + (y-nc).^2 );

    % Reduce high freq components by a power law f^-beta
    %to make a fractal surface
    F = F .* freq .^ -beta;

    % Set the DC level (= mean of the elevations) to zero
    F(nc,nc) = 0;	 

    % Take the inverse FFT
    M = real(ifft2(ifftshift(F))); 

    % clip out a portion with the requested dimensions
    M = M(1:ny,1:nx); 

    % scale elevations to a specified total variance
    M = M*sqrt(variance)/std(M(:));

    % % alternatively, scale elevations to produce a given relief, such that range(M(:)) = relief
    % M = M*relief/range(M(:)); 

    % % Optional: detrend the surface by subtracting a least-squares plane
    % M = detrend(M);

    
end
