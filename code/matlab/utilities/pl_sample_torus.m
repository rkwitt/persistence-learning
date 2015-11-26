function S = pl_sample_torus(N,r,R)
% SAMPLE_TORUS Draw a random sample  from a torus (in 3D).
%
%   S = PL_SAMPLE_TORUS(N,r,R) draws a sample of N points from a torus
%   defined by the inner radius r and the outer radius R.
%
% Author(s) Roland Kwitt, 2015

    phi = unifrnd(0,2*pi,N,1);
    theta = zeros(N,1);
    for i=1:N
        theta(i) = sample_theta(r,R);
    end
    
    p1 = (R+r*cos(theta)).*cos(phi);
    p2 = (R+r*cos(theta)).*sin(phi);
    p3 = r*sin(theta);
    
    S = [p1(:) p2(:) p3(:)];
end

function theta = sample_theta(r,R)
    while 1
        xvec = unifrnd(0,2*pi);
        yvec = unifrnd(0,1/pi);
        fx=(1+(r/R)*cos(xvec))/(2*pi);
        if (yvec<fx)
            break;
        end
    end 
    theta = xvec;
end