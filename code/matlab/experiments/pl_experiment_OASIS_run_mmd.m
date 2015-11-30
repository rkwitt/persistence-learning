function [K, p_values] = pl_experiment_OASIS_run_mmd(...
    mesh_mat_file, what, label, group_file, visit_file, options)
% PL_EXPERIMENT_OASIS_RUN_MMD runs MMD two-sample test
%
%   PL_EXPERIMENT_OASIS_RUN_MMD(MAT_FILE, WHAT, LABEL, GROUP_FILE, VISIT_FILE, OPTIONS) 
%   takes as input a MAT_FILE that is a saved struct with the following fields:
%
%       .config
%       .data
%
%       .config - Struct with fields
%
%           .alpha          - Scalar
%           .T1             - Array of HKS times
%           .target_scaling - Scaling of the mesh
%
%       .data - N x 1 cell array where each entry is a struct with the
%               following fields:
%
%           .V     - 3 x M matrix of vertex coordinates
%           .X     - M x 1 array of x-coordinates for each of the M vertices
%           .Y     - M x 1 array of y-coordinates for each of the M vertices 
%           .Z     - M x 1 array of z-coordinates for each of the M vertices
%           .TRIV  - K x 3 matrix of vertex indices of each mesh triangle 
%           .f_hks - M x len(T1) matrix of HKS times for each vertex
%      
%   Typically, MAT_FILE is the file that was saved when pre-processing the 
%   segmentations of the OASIS data. This is the same input file that is
%   used for PL_EXPERIMENT_OASIS_RUN_DIPHA
%
%   WHAT is a string which identifies the field to be loaded from the
%   MAT_FILE. LABEL is the prefix that is used for all output files.
%   GROUP_FILE and VISIT_FILE are files with meta-data information for all
%   subjects.
%
%   OPTIONS is a struct that is used to configure the MMD test. It needs to
%   have the following fields:
%
%       .dim                - Consider .dim homology, e.g., 1
%       .scales             - PSS kernel scales 10^-s, e.g., [1 2 3]
%       .trials             - Trials for bootstrapping, e.g., 10000
%       .alpha              - Significance level, e.g., 0.05
%       .collect_list_files - 0/1 
%       .compute_kernel     - 0/1
%       .run_mmd            - 0/1
%       .src_dir            - Source directory where .diagram files reside
%       .dst_dir            - Destination directory
%
% Author(s): Roland Kwitt, 2015

K =         -1; % Kernel matrices       
p_values =  -1; % MMD p-values

%--------------------------------------------------------------------------
%                                                                 Configure
%--------------------------------------------------------------------------
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
diagram_distance_binary = fullfile(root,'dipha-pss/build/diagram_distance');

groups = load(group_file);
visits = load(visit_file);
assert(length(groups)==length(visits));

% Demented vs. Non-Demented (on 1st visit)   
labels = zeros(length(groups),1);
p0 = find(groups==0 & visits==0); % 1st visit, non-demented subjects
p1 = find(groups==2 & visits==0); % 1st visit, demented subjects
labels(p0)=+1;
labels(p1)=-1;
pos = [p0(:); p1(:)];

% Load HKS data
tmp = load(mesh_mat_file); %#ok<NASGU>
eval(sprintf('X=tmp.%s', what));

% Create output directory
[s, ~, ~] = mkdir(options.dst_dir);
assert(s == 1);

%--------------------------------------------------------------------------
%                                                          Build list files
%--------------------------------------------------------------------------
if options.collect_list_files
    % Iterate over HKS times
    for t=1:length(X.config.T1)
        time_str = ['time_' num2str(t, '%.3d')];

        % File containing list of .diagram file names
        lst_file = fullfile(options.dst_dir, ...
            [...
                label '_' ...
                time_str '_lst.txt'
            ]); 

        % File containing labels of .diagram files
        lab_file = fullfile(options.dst_dir, ...
            [...
                label '_' ...
                time_str '_lab.txt'
            ]); 

        lst_fid = fopen(lst_file, 'w');
        lab_fid = fopen(lab_file, 'w');

        assert(lab_fid > 0);
        assert(lst_fid > 0);

        for i=1:length(pos)
            if isempty(X.data{pos(i)})
                continue;
            end
            [~,base_file_name,~] = fileparts(X.data{pos(i)}.file);
            diagram_file_name = fullfile(options.src_dir, ...
                [...
                    label '_' ...
                    base_file_name '_' ...
                    num2str(t, '%.3d') ...
                    '.diagram'
                ]);

            % Make sure file exists
            stat = exist(diagram_file_name, 'file' );
            assert(stat == 2);

            fprintf(lst_fid, '%s\n', diagram_file_name);
            fprintf(lab_fid, '%d\n', labels(pos(i)));
        end
        fclose(lst_fid);
        fclose(lab_fid);
    end
end

%--------------------------------------------------------------------------
%                                                        Kernel computation
%--------------------------------------------------------------------------
if options.compute_kernel
    % Iterate over HKS times
    for t=1:length(X.config.T1)
        time_str = ['time_' num2str(t, '%.3d')];
        lst_file = fullfile(options.dst_dir, ...
            [...
                label '_' ...
                time_str '_' ...
                'lst.txt'
             ]);

        % Iterate over PSS kernel time (sigma)
        for s=1:length(options.scales)
            fprintf('HKS-time [%d]: %.2f, PSS-sigma=%.2f\n', ...
                t, X.config.T1(t), options.scales(s));

            scale_str = ['scale_' num2str(s, '1e-%.3d')];
            dim_str = ['dim_' num2str(options.dim)]; 
            kernel_file = fullfile(options.dst_dir, ...
                [...
                    label '_' ...                   % e.g., cc
                    'K_inner_product' '_' ...       % IP kernel
                    time_str '_' ...                % HKS time
                    scale_str '_' ...               % PSS time (sigma)
                    dim_str '.txt'                  % Hom.-dim.
                ]);

            diagram_distance_options = [...
                ' --inner_product' ...
                ' --time ' num2str(10^-s,'%e') ...
                ' --dim ' num2str(options.dim) ' '];
            exec = [...
                diagram_distance_binary ...
                diagram_distance_options ...
                lst_file ...
                ' > ' ...
                kernel_file];
            system(exec);
        end
    end
end

%--------------------------------------------------------------------------
%                                                   Run MMD two-sample test
%--------------------------------------------------------------------------
if options.run_mmd
    K = cell(length(X.config.T1),1);
    dim_str = ['dim_' num2str(options.dim)];     

    % Iterate over HKS times
    for t=1:length(X.config.T1)
        time_str = ['time_' num2str(t, '%.3d')];

        % Iterate over PSS scales
        for s=1:length(options.scales)
            scale_str = ['scale_' num2str(s, '1e-%.3d')];
           
            % Construct kernel file name
            kernel_file = fullfile(options.dst_dir, ...
                [...
                    label '_' ...
                    'K_inner_product' '_' ...
                    time_str '_' ...
                    scale_str '_' ...
                    dim_str '.txt'
                ]);
            stat = exist(kernel_file, 'file' );
            assert(stat == 2);
            tmp = load(kernel_file);
            if s == 1
                K{t} = zeros(size(tmp,1), size(tmp,2), length(options.scales));
            end
            K{t}(:,:,s) = pl_normalize_kernel(tmp);
            K{t}(:,:,s) = exp(K{t}(:,:,s));
            
            fprintf('DONE preparing kernel (HKT time/PSS time)=(%.3f/%3.f)\n', ...
                X.config.T1(t), options.scales(s));
        end
    end

    % Diagnostics
    for t=1:length(X.config.T1)
        for s=1:length(options.scales)
            assert(sum(sum(isnan(K{t}(:,:,s))))==0);
        end
    end

    % Iterate over HKS times
    p_values = zeros(length(X.config.T1),length(options.scales));
    for t=1:length(X.config.T1)
        time_str = ['time_' num2str(t, '%.3d')];

        % Iterate over PSS scales
        for s=1:length(options.scales)
            KS = K{t}(:,:,s);

            % Construct label file name
            lab_file = fullfile(options.dst_dir, ...
                [...
                    label '_' ...
                    time_str '_lab.txt'
                ]); 
            
            y = load(lab_file);
            [~,idx] = sort(y);
            m = length(find(y==+1));
            n = length(find(y==-1));
            assert(n+m == length(y));
            KS = KS(idx, idx);

            % Call MMD 
            [testStat, ~, MMDarr] = pl_mmd(KS, m, n, ...
                options.trials, options.alpha);

            % Compute p-value
            p_values(t,s) = sum(MMDarr>testStat)/length(MMDarr);
            fprintf('DONE computing MMD (HKS time/PSS time)=(%.3f/%.3f)\n', ...
                X.config.T1(t), options.scales(s));
        end
    end
end