function pl_plot_pss( points, sigma, spacing, scaling )
% PL_PLOT_PSS plots the PSS mapping.
%
%   PL_PLOT_PSS(POINTS, SIGMA, SPACING, SCALING) plots the PSS map in the
%   form of a heat map for the PSS time parameter SIGMA. The grid to
%   evaluate the PSS map is provided via SPACING (i.e., the spacing of the
%   2D grid vertices in both birth and death direction). SCALING scales the
%   PSS map (typically set to 1 if you only have data from one diagram).
%
% Author(s): Roland Kwitt, 2015
    
    p = [points(:,1) points(:,2)];
    q = [points(:,2) points(:,1)];
        
    [X,Y] = meshgrid( spacing, spacing );
    gr = [X(:) Y(:)];
    f = zeros( size( gr, 1), 1 );
    for i=1:size(gr,1)
        f(i) = 1/scaling * eval_psi(gr(i,:), p, q, sigma);
    end
    
    %idx = gr(:,1)>=gr(:,2);
    %f(idx) = -Inf;
    
    m = size( X,1 );
    n = size( Y,1 );
    Z = reshape( f,[m n] );
    
    contourf(X,Y,Z,10,'LineWidth',1);
    axis off;
    axis equal;
    axis([-0.5 1.5 -0.5 1.5]);
    %mesh(X,Y,Z);
    %grid off;
    %box off;
    %axis off;
    %axis vis3d;
    %zoom(1.5);
end

function f = eval_psi( x, p, q, sigma )
    r = repmat(x, size( p,1 ), 1) - p;
    s = repmat(x, size( q,1 ), 1) - q;
    f = sum(exp(-sum( r.^2, 2 ) / (4*sigma)) - ...
    exp(-sum( s.^2, 2 ) / (4*sigma)));
end

