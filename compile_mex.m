% compile_mex.m
%
% Compiles Tadpole's C MEX files. Run it once from the Tadpole folder after
% downloading Tadpole, and again after changing any .c or .h file. It needs a C
% compiler that MATLAB supports:
%   Windows: Microsoft Visual C++ (free with Visual Studio Community) or the
%            MinGW-w64 add-on (MATLAB Add-On Explorer)
%   macOS:   Xcode or its command line tools (xcode-select --install)
%   Linux:   gcc
% Run "mex -setup C" to see or choose the compiler. If compiling fails with
% "cannot open file ... for writing", a MEX file is loaded in MATLAB: run
% "clear mex" and try again.
%
% The .c files include the header files (.h) in this folder, so compile from
% here.

files = {'mexD8', 'mexDinf', 'mexDms', ... % flow routing with depression filling
         'mexSPACE', ...                   % implicit river incision and sediment transport
         'mexHillPicard', ...              % implicit nonlinear hillslope diffusion
         'mexHillLinear'};                 % implicit linear hillslope diffusion

here = fileparts(mfilename('fullpath'));
olddir = cd(here);
cleanup = onCleanup(@() cd(olddir));
clear mex

for k = 1:numel(files)
    fprintf('Compiling %s.c ...\n', files{k});
    mex('-O', [files{k} '.c']);
end
fprintf('All MEX files compiled (%s).\n', mexext);
