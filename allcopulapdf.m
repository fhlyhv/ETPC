function PDF = allcopulapdf(typ,U,rho,nu)

% typ: 
%      elliptical
%      1 -- Gaussian;
%      2 -- t;
%      Archimedemean
%      3 -- Gumbel;
%      4 -- Clayton;
%      5 -- Frank;
%      Extreme value
%      6 -- Galambos;
%      7 -- HusterReiss;
%      8 -- tev

switch typ
    case 1 
        x = norminv(U);
        PDF = exp(-(rho.^2.*sum(x.^2,2)-2*rho.*prod(x,2))./(1-rho.^2)/2-log(1-rho.^2)/2);
    case 2
        x = tinv(U,repmat(nu,1,2));
        PDF = exp(gammaln((nu+2)/2)+gammaln(nu/2)-2*gammaln((nu+1)/2)-log(1-rho.^2)/2-(nu+2)/2.*...
            log(1+(sum(x.^2,2)-2*rho.*prod(x,2))./(nu.*(1-rho.^2)))+(nu+1)/2.*(log(1+x(:,1).^2./nu)+log(1+x(:,2).^2./nu)));
    case 3
        v = -log(U); % u is strictly in (0,1) => v strictly in (0,Inf)
        v = sort(v,2); vmin = v(:,1); vmax = v(:,2); % min/max, but avoid dropping NaNs
        nlogC = vmax.*(1+(vmin./vmax).^rho).^(1./rho);
        PDF = (rho - 1 + nlogC) ...
            .* exp(-nlogC + sum(repmat((rho-1),1,2).*log(v) + v, 2) + (1-2*rho).*log(nlogC));
        PDF(rho == 1) = 1;
    case 4
        logC = (-1./rho).*log(sum(U.^(-rho), 2) - 1);
        PDF = (rho+1) .* exp((2.*rho+1).*logC - sum((rho+1).*log(U), 2));
        PDF(rho==0) = 1;
    case 5
        expau = exp(rho .* U);
        PDF = -rho.*expm1(-rho) .* prod(expau,2) ...
            ./ (1 + exp(-rho).*prod(expau,2) - sum(expau, 2)).^2;
        PDF(rho==0) = 1;
    case 6
        PDF = Galamboscopulapdf(U,rho);
    case 7 
        PDF = HusterReisscopulapdf(U,rho);
    case 8 
        PDF = tevcopulapdf(U,rho,nu);
end