function P = pl_sample_linked_annuli( N, p, r1, r2, q, s1, s2, seed )
% PL_SAMPLE_LINKED_ANNULI draws a random sample of points from a double
% annulus (in 2D).
%
%   P = PL_SAMPLE_LINKED_ANNULI(N,P,R1,R2,Q,S1,S2,SEED) draws a random
%   sample of N 2D-points from the annuli defined by (P,R1,R2) and
%   (Q,S1,S2) where P and Q are the annuli centers and R1/S1 R2/S2 are the
%   inner and outer radii. SEED specifies the random number generator seed.
%
% Author(s): Roland Kwitt, 2015

A1 = pl_sample_annulus ( p, r1, r2, N, seed );
A2 = pl_sample_annulus ( q, s1, s2, N, seed );
P = [A1'; A2'];