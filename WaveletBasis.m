function [W,S,nscale,cscale] = WaveletBasis(Family,par,p)

% compute the wavelet basis matrix correponding to a specifict wavelet
% family
% Yu Hang, Jan, 2015, NTU

nscale = ceil(log2(p));
cscale = zeros(nscale,1);
W = [];
S = [];

for i = 1:nscale
    ns = nscale-i;
    cscale(i) = 2^ns;
    W1 = zeros(p,cscale(i));
    S1 = W1;
    for nc = 1:cscale(i)
        W1(:,nc) = MakeWavelet(ns,nc-1,Family,par,'Mother',p).';
        S1(:,nc) = MakeWavelet(ns,nc,Family,par,'Father',p).';
    end
    W1 = sparse(W1);
    S1 = sparse(S1);
    W = cat(2,W,W1);
    S = cat(2,S,S1);
end

cscale = [0;cumsum(cscale)];