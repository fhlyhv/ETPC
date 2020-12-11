function [rhohat, nuhat] = tevcopulafit(U)

[rhoi_mat, nui] = copulafit('t',U);
rhoi = rhoi_mat(1,2);
% nu_array = 1:round(nui)+5;
% nnu = length(nu_array);
% rho_array = zeros(1,nnu);
% nll_array = zeros(1,nnu);

lowerBndr = max(rhoi_mat(1,2)-0.1,0);
upperBndr = min(rhoi_mat(1,2)+0.1,1-1e-10);

lowerBndn = 0;
upperBndn = 10;
% opts = optimset(options);
nu_new = upperBndn;
rho_new = upperBndr;


while 1
    if abs(nu_new - nui)>1e-4 
        nu_new = fminbnd(@(nu)negloglike_tev(U(:,1),U(:,2),rhoi,nu),lowerBndn, upperBndn);
    end
    if  abs(rho_new - rhoi) <1e-6
        rho_new = fminbnd(@(rho)negloglike_tev(U(:,1),U(:,2),rho,nu_new),lowerBndr, upperBndr);
    end
    if abs(nu_new - nui)<1e-4 && abs(rho_new - rhoi) <1e-6
        break;
    else
        nui = nu_new;
        rhoi = rho_new;
    end
end

nuhat = nu_new;
rhohat = rho_new;

% for i = 1:nnu
%     [rho_array(i),nll_array(i)] = fminbnd(@(rho)negloglike_tev(U(:,1),U(:,2),rho,nu_array(i)),lowerBndr, upperBndr);
% end

% [~,id] = min(nll_array);
% rhohat = rho_array(id);
% nuhat = nu_array(id);


% ------------------------------------------------------------------------------
function nll = negloglike_tev(u1,u2,rho,nu)

logu1 = log(u1);
logu2 = log(u2);

c = gamma((nu+2)/2)/gamma((nu+1)/2)/sqrt((nu+1)*pi);
d = sqrt((nu+1)/(1-rho^2));

logq12 = logu1./logu2;
logq21 = logu2./logu1;

a = (logq12.^(1/nu)-rho)*d;
b = (logq21.^(1/nu)-rho)*d;

fa = 1+a.^2/(nu+1);
fb = 1+b.^2/(nu+1);

fl = logu1.*tcdf(a,nu+1)+logu2.*tcdf(b,nu+1);

logd = fl - logu1 - logu2 +...
    log((tcdf(a,nu+1)+d/nu*tpdf(a,nu+1).*logq12.^(1/nu)-d/nu*tpdf(b,nu+1).*logq21.^(1/nu+1)).*...
    (tcdf(b,nu+1)+d/nu*tpdf(b,nu+1).*logq21.^(1/nu)-d/nu*tpdf(a,nu+1).*logq12.^(1/nu+1))+...
    d/nu*((c*d*(nu+2)/nu/(nu+1)*a.*fa.^(-nu/2-2).*logq12.^(2/nu)-(1/nu+1)*tpdf(a,nu+1).*logq12.^(1/nu))./logu2+...
    (c*d*(nu+2)/nu/(nu+1)*b.*fb.^(-nu/2-2).*logq21.^(2/nu)-(1/nu+1)*tpdf(b,nu+1).*logq21.^(1/nu))./logu1));
    nll = -sum(logd);

