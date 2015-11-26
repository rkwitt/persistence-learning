function result = pl_experiment_pss_average(opt, out_dir)
% PL_EXPERIMENT_PSS_AVERAGE runs the PSS averaging as explained in the
% NIPS 2015 paper:
%
%   @inproceedings{Kwitt15a,
%       author    = {R.~Kwitt and S.~Huber and M.~Niethammer and W.~Lin and U.~Bauer},
%       title     = {Statistical Topological Data Analysis - A Kernel Perspective},
%       booktitle = {NIPS},
%       year      = {2015}}
%
% 
% Author(s): Roland Kwitt, Ulrich Bauer, 2015
        
if exist(out_dir,'dir')==0
    fprintf('Creating directory %s!\n', out_dir);
    mkdir(out_dir);
end

% DIPHA binary
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dipha_binary = fullfile(root,'dipha/build/dipha');

result.points = cell(opt.trials,1);
result.dmat =   cell(opt.trials,1);
result.pd =     cell(opt.trials,1); % Persistence diagrams

% Set RNG seed
rng(opt.seed);

% Draw sample from double annulus
FLA = pl_sample_linked_annuli( ...
    opt.N_points, ...
    opt.center1, ...
    opt.inner1, ...
    opt.outer1, ...
    opt.center2, ...
    opt.inner2, ...
    opt.outer2, ...
    -1 );

point_png = fullfile(out_dir, sprintf('FLA.png'));

% Plot all points
figure('visible','off');
plot(FLA(:,1),FLA(:,2),'r.','Markersize', 20 );
set(gcf, 'color', 'white');
set(gca, 'color', 'white');
axis tight;
axis off;
axis equal;
export_fig( point_png, '-r150' );    
close all;

all_pd_points = [];
for i=1:opt.trials
    fprintf( 'Trial %.3d/%.3d\n', i, opt.trials );
    
    % Sample points from full object
    sidx = randsample(2*10000, opt.N);
    
    % Get current sample
    points = FLA(sidx,:);
    
    % Plot current set of points
    point_png = fullfile(out_dir, sprintf('points_%.3d.png', i));
    figure('visible','off');
    plot(points(:,1),points(:,2),'r.','Markersize', 20 );
    set(gcf, 'color', 'white');
    set(gca, 'color', 'white');
    axis tight;
    axis off;
    axis equal;
    export_fig( point_png, '-r150' );    
    close all;
    
    % Pairwise Euclidean distance
    dmat = squareform( pdist( points ) );
    opt.max_val = 2*max( dmat(:) );
    dmat_file = fullfile( out_dir, sprintf( 'D_%.3d.bin', i ) );
    
    % Write distance matrix to DIPHA-compatible file
    save_distance_matrix(dmat, dmat_file);
    
    % Execute DIPHA using the distance matrix
    dipha_src_file = dmat_file;
    dipha_dst_file = fullfile( out_dir, sprintf( 'D_%.3d.pd', i ) );
    dipha_options = [...
        ' --upper_dim ' num2str( opt.max_dim ) ...
        sprintf(' %s', dipha_src_file ) ...
        sprintf(' %s', dipha_dst_file )];
    exec = ['/usr/local/bin/mpiexec -n 4 ' dipha_binary dipha_options ];
    system(exec);
    
    % Load persistence diagram by using DIPHA's MATLAB helpers
    [dim,b,d] = load_persistence_diagram( dipha_dst_file );
    data = [dim b d];
    data = data( dim==opt.target, 2:3 );
    txt_out = fullfile( out_dir, sprintf('pd_dim_%d_%.3d.txt', ...
        opt.target, i ) );
    dlmwrite( txt_out, data, 'delimiter', ' ' ); 
    
    result.pd{i} = [dim b d];
    result.points{i} = points;
    result.dmat{i} = dmat;
    
    % Iterate over a range of diagram_distance '--time' settings
    for s=1:length(opt.sigmas)
        png_file = fullfile(...
            out_dir, ...
            sprintf('pd_dim_%d_%.3d_sigma_%.5f.png', ...
                opt.target, ...
                i, ...
                opt.sigmas(s)));
        
        % Plot PSS mapping
        figure('visible','off');
        pl_plot_pss(data, opt.sigmas(s), -1.5:0.03:1.5, 1)
        set(gcf, 'color', 'white');
        set(gca, 'color', 'white');
        export_fig( png_file, '-r200' );
        close all;
    end
    all_pd_points = [all_pd_points; data];
end

for s=1:length(opt.sigmas)
    avg_png_file = fullfile( ...
        out_dir, ...
        sprintf('pd_dim_%d_avg_sigma_%.5f.png', ...
            opt.target, ...
            opt.sigmas(s)));
    figure('visible','off');
    pl_plot_pss(all_pd_points, opt.sigmas(s), -1.5:0.04:1.5, 1/opt.trials)
    set(gcf, 'color', 'white');
    set(gca, 'color', 'white');
    export_fig( avg_png_file, '-r200' );
    close all;
end
