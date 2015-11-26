function pl_write_persistence_diagram(filename,type,data)
    
    N = size(data,1);
    if strcmp(type,'dipha')
        filename = sprintf('%s.bin',filename);
        fid = fopen( filename, 'w' );
        fwrite( fid, 8067171840, 'int64' ); % DIPHA file
        fwrite( fid, 2, 'int64' ); % diagram ID
        fwrite( fid, N, 'int64' ); % pairs

        for i=1:size(data,1)
            fwrite( fid, data(i,1), 'int64');
            fwrite( fid, data(i,2), 'double');
            fwrite( fid, data(i,3), 'double');
        end
        fclose(fid);
    elseif strcmp(type ,'dionysus')
        filename = sprintf('%s.txt',filename);
        fid = fopen( filename, 'w' ); 
        p = find(data(:,1)==1);
        for i=1:length(p)
            fprintf(fid,'%.5f %.5f\n', data(p,2), data(p,3));
        end
        fclose(fid);
    else
        error('unknown type');
    end
end