function pl_setup()
% PL_SETUP Adds PL toolbox to the path.
%
% Authors: Roland Kwitt, 2015

[a,~,~] = fileparts(mfilename('fullpath')) ;
[a,~,~] = fileparts(a) ;
root = a ;

% These libraries we require for basic usage
addpath(fullfile(root, 'external/export_fig'        )) ;
addpath(fullfile(root, 'external/dipha/matlab'      )) ;
addpath(fullfile(root, 'matlab/experiments'         )) ;
addpath(fullfile(root, 'matlab/utilities'           )) ;
addpath(fullfile(root, 'matlab/testing'             )) ;

fprintf('Trying to load additional packages ...\n');

% The following libraries are only required for special purposes!
if (exist(fullfile(root,'external/sihks'), 'dir') == 7)
    addpath(fullfile(root, 'external/sihks'));
    fprintf('Found/Loaded: SIHKS\n');    
end
if (exist(fullfile(root,'external/STLRead'), 'dir') == 7)
    addpath(fullfile(root, 'external/STLRead'));
    fprintf('Found/Loaded: STLRead\n');        
end
if (exist(fullfile(root,'external/iso2mesh'), 'dir') == 7)
    addpath(fullfile(root, 'external/iso2mesh'));
    fprintf('Found/Loaded: iso2mesh\n');    
end

fprintf('PL ready.\n');