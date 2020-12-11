function PDF = Galamboscopulapdf(U,lambda)

u1 = U(:,1);
u2 = U(:,2);

nlogu1 = -log(u1);
nlogu2 = -log(u2);

fl = nlogu1.^(-1./lambda)+nlogu2.^(-1./lambda);

PDF = exp(fl.^(-lambda)).*...
    ((1-exp(-(lambda+1).*log(fl)-(1./lambda+1).*log(nlogu1))).*...
    (1-exp(-(lambda+1).*log(fl)-(1./lambda+1).*log(nlogu2)))+...
    (1+1./lambda).*exp(-(1./lambda+1).*log(nlogu1)-(lambda+2).*log(fl)-(1./lambda+1).*log(nlogu2)));