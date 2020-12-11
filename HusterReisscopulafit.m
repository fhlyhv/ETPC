function lambdahat = HusterReisscopulafit(U)

u1 = U(:,1);
u2 = U(:,2);

nloglf = @negloglike_HusterReiss;
lowerBnd = 0;
[lowerBnd,upperBnd] = bracket1D(nloglf,lowerBnd,5);

if ~isfinite(upperBnd)
    upperBnd = 1e100;
    %error(message('stats:copulafit:NoUpperBnd'));
end

lambdahat = fminbnd(nloglf, lowerBnd, upperBnd);


% ------------------------------------------------------------------------------
function nll = negloglike_HusterReiss(lambda)

logu1 = log(u1);
logu2 = log(u2);

logq12 = logu1./logu2;
logq21 = logu2./logu1;

a = lambda/2+log(logq12)/lambda;
b = lambda/2+log(logq21)/lambda;

fl = logu1.*normcdf(a)+logu2.*normcdf(b);

logd = fl-logu1-logu2+...
    log((normcdf(a)+normpdf(a)/lambda-normpdf(b).*logq21/lambda).*...
    (normcdf(b)+normpdf(b)/lambda-normpdf(a).*logq12/lambda)+...
    (normpdf(a)./logu2.*(a/lambda-1)+normpdf(b)./logu1.*(b/lambda-1))/lambda);

nll = -sum(logd);

end


% ------------------------------------------------------------------------------
function [nearBnd,farBnd] = bracket1D(nllFun,nearBnd,farStart)
% Bracket the minimizer of a (one-param) negative log-likelihood function.
% nearBnd is a point known to be a lower/upper bound for the minimizer,
% this will be updated to tighten the bound if possible.  farStart is the
% first trial point to test to see if it's an upper/lower bound for the
% minimizer.  farBnd will be the desired upper/lower bound.
bound = farStart;
upperLim = 1e12; % arbitrary finite limit for search
oldnll = nllFun(bound);
oldbound = bound;
while abs(bound) <= upperLim
    bound = 2*bound; % assumes lower start is < 0, upper is > 0
    nll = nllFun(bound);
    if nll > oldnll
        % The neg loglikelihood increased, we're on the far side of the
        % minimum, so the current point is the desired far bound.
        farBnd = bound;
        break;
    else
        % The neg loglikelihood continued to decrease, so the previous point
        % is on the near side of the minimum, update the near bound.
        nearBnd = oldbound;
    end
    oldnll = nll;
    oldbound = bound;
end
if abs(bound) > upperLim
    farBnd = NaN;
end
end

end