function P = ParabolicRidge(ny,nx,relief,beta,variance)

[X Y] = meshgrid(1:nx,1:ny);

P = -4*relief/(ny-1).^2 * (Y - 1 - (ny-1)/2).^2 + relief + RedNoise(ny,nx,beta,variance);
