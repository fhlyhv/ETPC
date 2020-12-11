% Jason K. Johnson
% MIT, June 2005
%
% SYNOPSIS: build an s x s GMRF based on "thin membrane" model which
% is useful to regularize image restoration problems by favoring
% images with small diffences between neighboring pixels.
%
% function [J,V] = thin_membrane(s1,s2)
% INPUTS dimensions s1 x s2 of 2D random field.
% OUTPUT let n = s1 * s2 (the number of vertices)
%   J: an n by n sparse spd information (inverse covariance) matrix.
%   V: an s1 by s2 full matrix containing indices of grid points in J.

%!!! column-major indexed instead of row-major

function [J,V,H] = thin_membrane(s1,s2) 

nv = s1*s2;

if nv < 2
    J = 1;
    V = 1;
    H = 1;
    return;
end;

D1 = ones(s1,3);
D1(:,2) = 0.0;
S1 = spdiags(D1,-1:1,s1,s1);

D2 = ones(s2,3);
D2(:,2) = 0.0;
S2 = spdiags(D2,-1:1,s2,s2);

I1 = speye(s1);
I2 = speye(s2);

% adjacency matrix (select storage format to minimize bandwidth of J)
%if (s1 <= s2)
  A = kron(I2,S1) + kron(S2,I1);
  V = reshape([1:nv],s1,s2);
% else
%   A = kron(I1,S2) + kron(S1,I2);
%   V = reshape([1:nv],s2,s1)';
% end

% thin membrane model
[ei,ej] = find(A);
kk = find(ei<ej);
ei = ei(kk);
ej = ej(kk);
ne = length(ei);
ee = [1:ne]';
H = sparse([ei;ej],[ee;ee],[ones(ne,1);-ones(ne,1)])';
J = H'*H;