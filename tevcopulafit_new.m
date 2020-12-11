%Trial tEV copula parameter estimation function

function [rhohat, nuhat] = tevcopulafit_new(U,option)

[rhoi_mat, nui] = multicopulafit('t',U);
rhoi = rhoi_mat(1,2);

lowerBndr = max(rhoi_mat(1,2)-0.1,0);
upperBndr = min(rhoi_mat(1,2)+0.1,1-1e-10);

lowerBndn = 0;
upperBndn = 30;

if option ==1
    rho_0 = rhoi;
    nu_0 =nui;
elseif option ==2
    rho_0 = rhoi;
    nu_0 = upperBndn;
elseif option ==3
    rho_0 = upperBndr;
    nu_0 =nui;
elseif option ==4
    rho_0 = upperBndr;
    nu_0 = upperBndn;
elseif option ==5
    rho_0 = rhoi;
    nu_0 = lowerBndn;
elseif option ==6
    rho_0 = upperBndr;
    nu_0 = lowerBndn;
end

options =  optimset('Display','off','Algorithm','interior-point');

answer = fmincon(@(x)negloglike_tev(x,U(:,1),U(:,2)),[rho_0;nu_0],[],[],[],[],[lowerBndr;lowerBndn],[upperBndr;upperBndn],[],options);

rhohat = answer(1);
nuhat = answer(2);

function nll = negloglike_tev(x,u1,u2)

logu1 = log(u1);
logu2 = log(u2);

c = gamma((x(2)+2)/2)/gamma((x(2)+1)/2)/sqrt((x(2)+1)*pi);
d = sqrt((x(2)+1)/(1-x(1)^2));

logq12 = logu1./logu2;
logq21 = logu2./logu1;

a = (logq12.^(1/x(2))-x(1))*d;
b = (logq21.^(1/x(2))-x(1))*d;

fa = 1+a.^2/(x(2)+1);
fb = 1+b.^2/(x(2)+1);

fl = logu1.*tcdf(a,x(2)+1)+logu2.*tcdf(b,x(2)+1);

logd = fl - logu1 - logu2 +...
    log((tcdf(a,x(2)+1)+d/x(2)*tpdf(a,x(2)+1).*logq12.^(1/x(2))-d/x(2)*tpdf(b,x(2)+1).*logq21.^(1/x(2)+1)).*...
    (tcdf(b,x(2)+1)+d/x(2)*tpdf(b,x(2)+1).*logq21.^(1/x(2))-d/x(2)*tpdf(a,x(2)+1).*logq12.^(1/x(2)+1))+...
    d/x(2)*((c*d*(x(2)+2)/x(2)/(x(2)+1)*a.*fa.^(-x(2)/2-2).*logq12.^(2/x(2))-(1/x(2)+1)*tpdf(a,x(2)+1).*logq12.^(1/x(2)))./logu2+...
    (c*d*(x(2)+2)/x(2)/(x(2)+1)*b.*fb.^(-x(2)/2-2).*logq21.^(2/x(2))-(1/x(2)+1)*tpdf(b,x(2)+1).*logq21.^(1/x(2)))./logu1));
    
nll = -sum(logd);