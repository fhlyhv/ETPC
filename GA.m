function B = GA(Jp,cplmat,W1,S1,W2,S2,nscale,cscale1,cscale2,a1,a2,tol)

% computing optimal edge weigth matrix B using stochastic gradient ascent algorithm
% Yu Hang, Jan. 2015, NTU

% initialize B mattrix
p = size(Jp,1);
n = size(cplmat,1);
[edgerow,edgecol] = find(triu(Jp,1));
nedge = length(edgerow);

W0 = log(a2*(ones(nedge,1))/sqrt(nedge));
W = W0;

% initialize parameters for ADADELTA
eps = 1e-6;
rho = 0.95;
d2W = 0;
g2W = 0;
gW_array = zeros(nedge,n);
% eta = 1e-3;

% compute matrix E
E = sparse(nedge,p);
for i = 1:nedge
    if edgecol(i)<p
        E(i,edgerow(i)) = 1;
        E(i,edgecol(i)) = -1;
    else
        E(i,edgerow(i)) = 1;
    end
end

% begin stochastic gradient ascent algorithm
for nt = 1:5e5
    
    Bv = exp(W);
%     ne = ceil(n*rand(1));
    
    % compute low rank approximation matrix
%     L = LowRankApp(4,W1,S1,W2,S2,nscale,cscale1,cscale2);
%     RC = QBC\L;
    

%     gW = (n*cplmat(ne,:).'.*sum((E*RC).*(E*L),2)-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
    parfor ne = 1:n
        BC = sparse(edgerow,edgecol,Bv.*cplmat(ne,:).',p,p);
        BC = BC+BC.';
        QBC = spdiags(sum(BC,2),0,-BC);
        QBC = [QBC(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];
        gW_array(:,ne) = cplmat(ne,:).'.*sum(E*(QBC\speye(p)).*E,2);
    end
    gW = (mean(gW_array,2)-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
    g2W = rho*g2W+(1-rho)*gW.^2;
    dW =  sqrt(d2W+eps)./sqrt(g2W+eps).*gW;
    d2W = rho*d2W+(1-rho)*dW.^2;
    W = W+dW;%eta*gW; %
    
    if rem(nt,10) == 0
        Wh = W;
        save resulttmptg2;
        mean(abs(Wh-W0))
        if mean(abs(Wh-W0)) < 1e-4
            break;
        else
            W0 = Wh;
%             eta = eta*0.995;
        end
    end
end


W0 = log(a2*exp(W)/norm(exp(W)));
W = W0;
d2W = 0;
g2W = 0;
% eta = 1e-3;

for nt = 1:5e5
    
    Bv = exp(W);
    B = sparse(edgerow,edgecol,Bv,p,p);
    B = B+B.';
    QB = spdiags(sum(B,2),0,-B);
    QB = [QB(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];
    % compute low rank approximation matrix
%     L = LowRankApp(4,W1,S1,W2,S2,nscale,cscale1,cscale2);
%     R = QB\L;
    
%     ne = ceil(n*rand(1));
    
%     RC = QBC\L;
    
%     gW = (n*sum((spdiags(cplmat(ne,:).',0,nedge,nedge)*(E*RC)-E*R).*(E*L),2)-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
    parfor ne = 1:n
        BC = sparse(edgerow,edgecol,Bv.*cplmat(ne,:).',p,p);
        BC = BC+BC.';
        QBC = spdiags(sum(BC,2),0,-BC);
        QBC = [QBC(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];
        gW_array(:,ne) = cplmat(ne,:).'.*sum(E*(QBC\speye(p)).*E,2);
    end
    gW = (mean(gW_array,2)-sum(E*(QB\speye(p)).*E,2)-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
    g2W = rho*g2W+(1-rho)*gW.^2;
    dW =  sqrt(d2W+eps)./sqrt(g2W+eps).*gW;
    d2W = rho*d2W+(1-rho)*dW.^2;
    W = W+dW; %eta*gW; %
    
    if rem(nt,10) == 0
        Wh = W;
        save resulttmptg2;
        mean(abs(Wh-W0))
        if mean(abs(Wh-W0)) < 1e-6
            break;
        else
            W0 = Wh;
%             eta = eta*0.995;
        end
    end
end

B = sparse(edgerow,edgecol,a2/norm(exp(Wh))*exp(Wh),p,p);
B = B+B.';
    
    
    
    