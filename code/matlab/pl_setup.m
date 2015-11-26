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

fprintf('PL ready.\n');

end