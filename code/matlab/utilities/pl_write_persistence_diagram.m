function pl_write_persistence_diagram(filename,type,data)
% PL_WRITE_PERSISTENCE_DIAGRAM writes points in a persistence diagram to a
% file compatible with DIPHA or DIONYSUS.
%
%   PL_WRITE_PERSISTENCE_DIAGRAM(FILENAME, TYPE, DATA) writes points of a
%   persistence diagram in the N x 3 array DATA to a file (FILENAME) that
%   can be read by either DIPHA (when TYPE is 'dipha') or DIONYSUS (when 
%   TYPE is 'dionysus'). 
%
%   The DATA array is contains in the FIRST column the dimension of the
%   points, e.g., 0, 1, the SECOND column contains the birth times, the
%   THIRD column contains the death times. For example,
%
%   data = [
%       0 1.1 1.2, ...
%       0 1.3 1.4];
%   pl_write_persistence_diagram('/tmp/diagram','dipha',data); 
%
%   will write the file /tmp/diagram.bin in DIPHA-compatible format
%   containing two points (1.1,1.2) and (1.3,1.4).
%   
% Author(s): Roland Kwitt, 2015

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