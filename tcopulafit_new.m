function [rhohat, nuhat] = tcopulafit_new(U)

% [rhomat,~] = copulafit('t',U);
% rhohat = rhomat(1,2);

[rhohat,~] = corr(U(:,1),U(:,2),'type','Kendall');
rhohat = sin(pi*rhohat/2);

lowerBndn = 0;
upperBndn = 10;

options =  optimset('Display','off','Algorithm','interior-point');

nu_0 = 100;

%nuhat = fmincon(@(nu)negloglike_t(rhohat,nu,U(:,1),U(:,2)),nu_0,[],[],[],[],lowerBndn,upperBndn,[],options);
nuhat = fminbnd(@(nu)negloglike_t(rhohat,nu,U(:,1),U(:,2)),lowerBndn,upperBndn,options);

function nll = negloglike_t(rho,nu,u1,u2)

const = gamma((nu+2)/2)/gamma(nu/2)/nu/sqrt(1-rho^2)/pi;

x1 = tinv(u1,nu);
x2 = tinv(u2,nu);

logd = log(const) - log(tpdf(x1,nu)) - log(tpdf(x2,nu)) - 0.5*(nu+1)*log(1 + (x1.^2 + x2.^2 - 2*rho*x1.*x2)/nu/(1-rho)^2);
    
nll = -sum(logd);