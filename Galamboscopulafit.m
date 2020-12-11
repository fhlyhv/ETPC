function lambdahat = Galamboscopulafit(U)

% lambda = 1/theta

u1 = U(:,1);
u2 = U(:,2);

nloglf = @negloglike_Gamlambos;
lowerBnd = 0;
[lowerBnd,upperBnd] = bracket1D(nloglf,lowerBnd,5);

if ~isfinite(upperBnd)
    error(message('stats:copulafit:NoUpperBnd'));
%     upperBnd = [];
end

lambdahat = real(fminbnd(nloglf, lowerBnd, upperBnd)); %

% --------------------------------------------------------------------
function nll = negloglike_Gamlambos(lambda)

nlogu1 = -log(u1);
nlogu2 = -log(u2);

fl = nlogu1.^(-1/lambda)+nlogu2.^(-1/lambda);

logd = fl.^(-lambda)+ ...
    log((1-exp(-(lambda+1)*log(fl)-(1/lambda+1)*log(nlogu1))).*...
    (1-exp(-(lambda+1)*log(fl)-(1/lambda+1)*log(nlogu2)))+...
    (1+1/lambda)*exp(-(1/lambda+1)*log(nlogu1)-(lambda+2)*log(fl)-(1/lambda+1)*log(nlogu2)));

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