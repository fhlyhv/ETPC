function PDF = tevcopulapdf(U,rho,nu)

u1 = U(:,1);
u2 = U(:,2);

logu1 = log(u1);
logu2 = log(u2);

c = gamma((nu+2)/2)./gamma((nu+1)/2)./sqrt((nu+1)*pi);
d = sqrt((nu+1)./(1-rho.^2));

logq12 = logu1./logu2;
logq21 = logu2./logu1;

a = (logq12.^(1./nu)-rho).*d;
b = (logq21.^(1./nu)-rho).*d;

fa = 1+a.^2./(nu+1);
fb = 1+b.^2./(nu+1);

fl = logu1.*tcdf(a,nu+1)+logu2.*tcdf(b,nu+1);

PDF = exp(fl)./u1./u2.*((tcdf(a,nu+1)+d./nu.*tpdf(a,nu+1).*logq12.^(1./nu)-d./nu.*tpdf(b,nu+1).*logq21.^(1./nu+1)).*...
    (tcdf(b,nu+1)+d./nu.*tpdf(b,nu+1).*logq21.^(1./nu)-d./nu.*tpdf(a,nu+1).*logq12.^(1./nu+1))+...
    d./nu.*((c.*d.*(nu+2)./nu./(nu+1).*a.*fa.^(-nu/2-2).*logq12.^(2./nu)-(1./nu+1).*tpdf(a,nu+1).*logq12.^(1./nu))./logu2+...
    (c.*d.*(nu+2)./nu./(nu+1).*b.*fb.^(-nu/2-2).*logq21.^(2./nu)-(1./nu+1).*tpdf(b,nu+1).*logq21.^(1./nu))./logu1));