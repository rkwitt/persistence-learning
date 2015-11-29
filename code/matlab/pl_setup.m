function pl_setup()
% PL_SETUP Adds PL toolbox to the path.
%
% Authors: Roland Kwitt, 2015

[a,~,~] = fileparts(mfilename('fullpath')) ;
[a,~,~] = fileparts(a) ;
root = a ;

addpath(fullfile(root, 'export_fig'        )) ;
addpath(fullfile(root, 'dipha/matlab'      )) ;
addpath(fullfile(root, 'matlab/experiments')) ;
addpath(fullfile(root, 'matlab/utilities'  )) ;
addpath(fullfile(root, 'matlab/testing'    )) ;

if (exist(fullfile(root,'external/sihks'), 'dir') == 7)
    addpath(fullfile(root, 'external/sihks'));
end
if (exist(fullfile(root,'external/STLRead'), 'dir') == 7)
    addpath(fullfile(root, 'external/STLRead'));
end
if (exist(fullfile(root,'external/iso2mesh'), 'dir') == 7)
    addpath(fullfile(root, 'external/iso2mesh'));
end
fprintf('PL ready.\n');
end