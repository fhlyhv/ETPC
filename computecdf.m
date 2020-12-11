function Fx = computecdf(XDat,U,G,S,q)

% compute cdf all the data
% for those data above the threshold, use GP
% for those data under the threshold, use empirical cdf

[n,p] = size(XDat);
Fx = zeros(n,p);

parfor i = 1:p
    X = XDat(:,i);
    f = zeros(n,1);
    f(X>=U(i)) = q+(1-q)*gpcdf(X(X>=U(i)),G(i),S(i),U(i));
    f(X<U(i)) = q*empcdf_con(X(X<U(i)));
    Fx(:,i) = f;
end
