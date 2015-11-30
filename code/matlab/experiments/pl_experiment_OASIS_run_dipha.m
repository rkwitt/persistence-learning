function pl_experiment_OASIS_run_dipha(mat_file, what, label, out_dir)
% PL_EXPERIMENT_RUN_DIPHA constructs simplicial complexes from meshes with
% Heat-Kernel signature computed at each vertex of the mesh.
%
%   PL_EXPERIMENT_OASIS_RUN_DIPHA(MAT_FILE, WHAT, LABEL, OUT_DIR) takes as
%   input a MAT_FILE that is a saved struct with the following fields:
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
%   segmentations of the OASIS data.
%
%   WHAT is a string which identifies the field to be loaded from the
%   MAT_FILE. LABEL is the prefix that is used for all output files and
%   OUT_DIR is the directory where all output files will be written to
%   (will be created if it does not exist).
%
% Author(s): Roland Kwitt

    % Load data
    Y = load(mat_file);
    eval(sprintf('Y=Y.%s', what));
    
    % DIPHA binary
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    dipha_binary = fullfile(root,'external/dipha/build/dipha');
    [status,~,~] = mkdir(out_dir);
    assert(status == 1);
    
    %----------------------------------------------------------------------
    %                                   Create complexes + compute diagrams
    %----------------------------------------------------------------------
    for i=1:length(Y.data)
        X = Y.data{i};
        num_times = length(Y.config.T1);
        
        if isempty(X)
            continue;
        end
       
        % Iterate over HKS times
        for t=1:num_times
            triangles = X.TRIV;
            vertex_values = X.f_hks(:,t);
            [~,base_file_name,~] = fileparts(X.file);

            % Write simplicial complexes
            complex_file_name = fullfile(out_dir, ...
                [...
                    label '_' ...
                    base_file_name '_' ...
                    num2str(t, '%.3d') '.complex'
                ]);

            % Avoid unnecessary recomputation
            if (exist(complex_file_name, 'file' ) ~= 2)
                pl_mesh2dipha( ...
                    triangles, ...
                    vertex_values, ...
                    complex_file_name);
            end

            % Compute persistence diagrams
            diagram_file_name = fullfile(out_dir, ...
                [...
                    label '_' ...
                    base_file_name '_' ...
                    num2str(t, '%.3d') '.diagram' ...
                ]);
            exec = [dipha_binary ' ' ...
                complex_file_name ' ' ...
                diagram_file_name];

            % Again, avoid unnecessary recomputation
            if (exist(diagram_file_name, 'file' ) ~= 2)
                system( exec );
            end
        end
        
        fprintf('DONE with %s\n', X.file);
    end
end
