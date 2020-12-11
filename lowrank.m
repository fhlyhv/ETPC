function B = lowrank(Family,par,p1,p2,r,nscale)

% generate low rank matrix B for 2D problem
% Family -- wavelet family  p1 -- no. of rows  p2 -- no of columns  r --
% rank of each scale

% Yu Hang, NTU, Jan, 2015


nscale1 = ceil(log2(p1));
nscale2 = ceil(log2(p2));
if nargin <6
    nscale = min(nscale1,nscale2);  % match from the finest scale
end

B = [];

for ns = 1:nscale
    ns1 = nscale1-ns;
    pc1 = 2^ns1;
    W1 = zeros(p1,pc1);
    S1 = W1;
    for nc = 1:pc1
        W1(:,nc) = MakeWavelet(ns1,nc-1,Family,par,'Mother',p1).';
        S1(:,nc) = MakeWavelet(ns1,nc,Family,par,'Father',p1).';
    end
    W1 = sparse(W1);
    S1 = sparse(S1);
    if pc1 > r
        C = spdiags(sign(binornd(1,0.5,r,pc1/r)-0.5),0:-r:1-pc1,pc1,r);
        W1 = W1*C;
        S1 = S1*C;
    end
    
    
    ns2 = nscale2-ns;
    pc2 = 2^ns2;
    W2 = zeros(p2,pc2);
    S2 = W2;
    for nc = 1:pc2
        W2(:,nc) = MakeWavelet(ns2,nc-1,'Coiflet',par,'Mother',p2).';
        S2(:,nc) = MakeWavelet(ns2,nc,'Coiflet',par,'Father',p2).';
    end
    W2 = sparse(W2);
    S2 = sparse(S2);
    if pc2 > r
        C = spdiags(sign(binornd(1,0.5,r,pc2/r)-0.5),0:-r:1-pc2,pc2,r);
        W2 = W2*C;
        S2 = S2*C;
    end
    
    B = cat(2,B,[kron(W2,W1),kron(W2,S1),kron(S2,W1)]);
    
end
B = cat(2,B,kron(S2,S1));
    


