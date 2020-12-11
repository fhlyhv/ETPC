function SS=cov_normalize(S)

A = sqrt(diag(S));
SS = S./(A*A.');
    