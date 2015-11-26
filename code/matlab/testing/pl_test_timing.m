function stat = pl_test_timing(out_dir, range, M, dim, time)
% PL_TEST_TIMING Runs timing tests for kernel with/without the
% Fast-Gauss-Transform.
%
%   STAT = PL_TEST_TIMING(OUT_DIR, RANGE, M, DIM, TIME) runs kernel timing
%   tests for the PSS kernel.
%
%   For each value K in RANGE, the function creates M persistence diagrams
%   with K random points and writes them into OUT_DIR which is created
%   in case it does not exist. It then runs the kernel (1) without the
%   Fast-Gauss-Transform (FGT) and (2) with the Fast-Gauss-Transform using 
%   DIM and TIME as kernel parameters (see ./diagram_distance --help for 
%   details).
%   
%   The output in STAT is a length(RANGE) x 3 array where each row 
%   corresponds to the number of points in the diagram. The FIRST column
%   is the execution time (in seconds) for the kernel without FGT. 
%   The SECOND column is the executing time (in seconds) for the kernel
%   with the Fast-Gauss-Transform and the THIRD column is the Frobenius
%   norm between the generated kernel (Gram) matrices.
%
% Author(s): Roland Kwitt, 2015

    [~,~,~] = mkdir(out_dir);
    
    % Binary
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    diagram_distance = fullfile(root,'dipha-pss/build/diagram_distance');
    
    % Time storage 
    % idx=1: WITHOUT FGT
    % idx=2: WITH FGT
    % idx=3: Difference in Frobenius norm
    stat = zeros(length(range),3);
    
    for i=1:length(range)
        N = range(i);
        fprintf('Testing for %.3d points\n', N);
        
        % File to hold diagram file names
        list_file = fullfile(out_dir, sprintf('list_%.3d.txt', N));
        fid = fopen(list_file, 'w');
        
        % Number of diagrams to create
        for j=1:M
            b = rand(N,1);      % birth time
            d = b + rand(N,1);  % death time
            data = [ones(N,1)*dim b d];
            
            out_file = sprintf('test_%.3d_%.3d', N, j);
            out_file = fullfile(out_dir, out_file);
            pl_write_persistence_diagram(out_file, 'dipha', data);
            
            fprintf(fid, '%s.bin\n', out_file);
        end
        
        fclose(fid);
        
        % GRAM matrices for comparison
        gram_matrix_file_nIFGT = fullfile(out_dir, ...
            sprintf('K_nIFGT_%.3d.txt', N));
        gram_matrix_file_wIFGT = fullfile(out_dir, ...
            sprintf('K_wIFGT_%.3d.txt', N));
        
        % (1) Without FGT
        options = ['--inner_product --time ' ...
            num2str(time,'%e') ...
            ' --dim ' num2str(dim) ' '];
        exec = [diagram_distance ' ' options list_file ' > ' gram_matrix_file_nIFGT];
        tic; system( exec ); ela_nIFGT = toc;
        
        % (2) With FGT
        options = ['--inner_product --time ' ...
            num2str(time,'%e') ...
            ' --use_fgt ' ...
            ' --dim ' num2str(dim) ' '];
        exec = [diagram_distance ' ' options list_file ' > ' gram_matrix_file_wIFGT];
        tic; system( exec ); ela_wIFGT = toc;
        
        K1 = load(gram_matrix_file_nIFGT); % GRAM matrix without FGT
        K2 = load(gram_matrix_file_wIFGT); % GRAM matrix with FGT
        
        stat(i,1) = ela_nIFGT;
        stat(i,2) = ela_wIFGT;
        stat(i,3) = norm(K1-K2,'fro');
    end
end