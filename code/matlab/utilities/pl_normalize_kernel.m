% normalize_kernel normalizes a kernel Gram matrix
%
% This is standard postprocessing, described in the paper 'Scalable kernels 
% for graphs with continuous attributes (A. Feragen, N. Kasenburg, 
% J. Petersen, M. de Bruijne, K. Borgwardt), Neural Information Processing 
% Systems (NIPS) 2013'.
%
% K_norm = normalize_kernel(K)
%
% Computes the normalized Gram matrix K_norm from an unnormalized 
% Gram matrix K by setting K_norm(i,j) = K(i,j)/sqrt(K(i,i)K(j,j)). 
%
% Now all diagonal entries in K_norm are 1.

% Copyright (C) 2013 Aasa Feragen
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
% Author: Aasa Feragen <aasa at diku dot dk>
% 
% 2013-10-26 Aasa Feragen <aasa at diku dot dk>
% * Initial revision


function K_norm = pl_normalize_kernel(K)   
    Kdiag = diag(K);
    normalizers = sqrt(Kdiag*Kdiag');
    K_norm = K./normalizers;    
end
    