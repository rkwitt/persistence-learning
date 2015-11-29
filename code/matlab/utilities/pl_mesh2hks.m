function result = pl_mesh2hks(base_dir, subjects, target_scaling, varargin)
% PL_MESH2HKS computes heat-kernel signatures (HKS) for 3D meshes of 
% available in STL format.
%
% Author(s): Roland Kwitt, 2015

    % Additional arguments
    if nargin == 3
       alpha = 2;
       T1 = 1:0.5:10.5;
       verbose = 1;
    elseif nargin == 4
       alpha = varargin{1};
       T1 = 1:0.5:10.5;
       verbose = 1;
    elseif nargin == 5
        alpha =  varargin{1};
        T1 = varargin{2};
        verbose = 1;
    end

    % Create cell array for data from all subjects
    subject_data = cell(length(subjects),1);

    % Iterate over all subjects
    for subject_id=1:length(subjects)
        file_name = subjects{subject_id};
        
        % Sanity check:
        %
        % Make sure mesh file exists; if not, it might have been removed
        stat = exist(fullfile(base_dir, file_name), 'file' );
        if stat ~= 2
            if verbose
                fprintf('Subject %d missing or removed\n', subject_id);
            end
            continue;
        end

        % Possibly repair mesh
        [faces, vertices] = stlread(fullfile(base_dir, file_name));
        [vertices_fixed, faces_fixed] = ...
            meshcheckrepair(vertices, faces);

        % Normalize meshes
        V = vertices_fixed;
        V = V - repmat(mean(V),size(V,1),1);
        V = V/norm(V);
        V = target_scaling*V;

        % Save (scaled) mesh data
        subject_data{subject_id}.V = V';
        subject_data{subject_id}.X = V(:,1)';
        subject_data{subject_id}.Y = V(:,2)';
        subject_data{subject_id}.Z = V(:,3)';
        subject_data{subject_id}.TRIV = faces_fixed;
        subject_data{subject_id}.file = fullfile(base_dir, file_name);

        % Compute HKS
        subject_data{subject_id}.f_hks = compute_hks(...
            subject_data{subject_id}, alpha, T1);
    end

    result.config.target_scaling = target_scaling;
    result.config.alpha = alpha;
    result.config.T1 = T1;
    result.data = subject_data;
end

function f_hks = compute_hks(data,alpha,T1)
    opt.dtype = 'cotangent';
    [W,A] = mshlp_matrix(data, opt);
    A = spdiags(A, 0, size(A,1), size(A,1));
    [evecs,evals] = eigs(W, A, 50, 'SM');
    evals = -diag(evals);
    f_hks = hks(evecs, evals, alpha.^T1);
end
