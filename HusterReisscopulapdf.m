function PDF = HusterReisscopulapdf(U,lambda)

u1 = U(:,1);
u2 = U(:,2);



logu1 = log(u1);
logu2 = log(u2);

logq12 = logu1./logu2;
logq21 = logu2./logu1;

a = lambda/2+log(logq12)./lambda;
b = lambda/2+log(logq21)./lambda;

% if ~isreal(a) || ~isreal(b)
%     a = real(a);
%     b = real(b);
% end

fl = logu1.*normcdf(a)+logu2.*normcdf(b);

PDF = exp(fl)./u1./u2.*((normcdf(a)+normpdf(a)./lambda-normpdf(b).*logq21./lambda).*...
    (normcdf(b)+normpdf(b)./lambda-normpdf(a).*logq12./lambda)+...
    (normpdf(a)./logu2.*(a./lambda-1)+normpdf(b)./logu1.*(b./lambda-1))./lambda);