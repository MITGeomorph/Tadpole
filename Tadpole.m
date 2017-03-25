function solution = Tadpole(initial,p)

% Tadpole.m
%
% solution = Tadpole(initial,p)
%

%   Output argument:
%
%   solution     The elevation matrix after the final timestep
%
%
%   Input arguments:
%
%   initial      Matrix of initial elevations
%
%   p            Struct array of parameters that includes:
%
%  Two struct arrays are used to handle all data: 
%   p contains parameters
%   g contains grids
%     

% check arguments

if nargin ~= 2, help adm; return; end


% run model

[p g] = TadpoleInitialize(initial,p);


[p g] = TadpoleRun(p,g);


[p g] = TadpoleFinalize(p,g);

            
% assign solution to output argument
solution = g.output;