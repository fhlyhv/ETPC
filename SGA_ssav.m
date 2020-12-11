function B = SGA_ssav(Jp,cplmat,a1,a2)

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

gW_array = zeros(nedge,n);

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

c = graphcolor_irregular(Jp^4);

% quasi random index
ida= haltonset(1,'Skip',1e3,'Leap',1e2);
ida = scramble(ida,'RR2'); 
ida = qrandstream(ida);

% begin stochastic gradient ascent algorithm
flag = 1;
nid = 1:n;
g = 0;
h = 0;
tau = n;
eta = 0.05;
rho = 0.01;
eps = 1e-6;
dW = 0;
dW1 = 0;

for nt = 1:30*n  
    
    Bv = exp(W);
    ne = ceil(n*qrand(ida,1));
    BC = sparse(edgerow,edgecol,Bv.*cplmat(ne,:).',p,p);
    BC = BC+BC.';
    QBC = spdiags(sum(BC,2),0,-BC);
    QBC = [QBC(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];
    % compute low rank approximation matrix
    L = sparse((1:p).',c,sign(randi([0,1],p,1)-0.5)); %LowRankApp(4,W1,S1,W2,S2,nscale,cscale1,cscale2);
    RC = QBC\L;
    
    
    gW_array(:,ne) = cplmat(ne,:).'.*sum((E*RC).*(E*L),2); %cplmat(ne,:).'.*sum(E*(QBC\speye(p)).*E,2);
    
    if flag == 0
        gW = (mean(gW_array,2)-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
    else
        nid = setdiff(nid,ne);
        gW = (sum(gW_array,2)/(n-length(nid))-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
        if isempty(nid)
            flag = 0;
        end
    end
    
%     if nt == 1
%         dW2 = gW.^2;
%     else
%         dW2 = (1-eta)*dW2 + eta*gW.^2;
%     end

    W = W+rho*gW;
    
    
    if rem(nt,3*n) == 0
        Wh = W;
        mean(abs(Wh-W0))
        if mean(abs(Wh-W0)) < 1e-3
            break;
        else
            W0 = Wh;
            rho = rho*0.6;
        end
    end
end


W0 = log(a2*exp(W)/norm(exp(W)));
W = W0;
gW_array = zeros(nedge,n);
nid = 1:n;
flag = 1;
g = 0;
h = 0;
% tau = n;
Wa = [];
% dW = 0;
% dW1 = 0;
tau = n;
rho = 0.01;

for nt = 1:5e5
    
    Bv = exp(W);
%     if norm(Bv) < 8
%         norm(Bv)
%     end
    B = sparse(edgerow,edgecol,Bv,p,p);
    B = B+B.';
    QB = spdiags(sum(B,2),0,-B);
    QB = [QB(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];

    % compute low rank approximation matrix
    L = sparse((1:p).',c,sign(randi([0,1],p,1)-0.5)); %LowRankApp(4,W1,S1,W2,S2,nscale,cscale1,cscale2);
    R = QB\L;
    
    ne = ceil(n*qrand(ida,1));
    BC = sparse(edgerow,edgecol,Bv.*cplmat(ne,:).',p,p);
    BC = BC+BC.';
    QBC = spdiags(sum(BC,2),0,-BC);
    QBC = [QBC(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];
    RC = QBC\L;
    
    gW_array(:,ne) = sum((spdiags(cplmat(ne,:).',0,nedge,nedge)*(E*RC)-E*R).*(E*L),2);
    %cplmat(ne,:).'.*sum(E*(QBC\speye(p)).*E,2)-sum(E*(QB\speye(p)).*E,2);
    if flag == 0
        gW = (mean(gW_array,2)-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
    else
        nid = setdiff(nid,ne);
        gW = (sum(gW_array,2)/(n-length(nid))-4*a1*(norm(Bv)^2-a2^2)*Bv).*Bv;
        if isempty(nid)
            flag = 0;
        end
    end
    
%     if nt == 1
%         dW2 = gW.^2;
%     else
%         dW2 = (1-eta)*dW2 + eta*gW.^2;
%     end

    W = W+rho*gW;  %./sqrt(dW2+eps)
    Wa = cat(2,Wa,W);
    
    if rem(nt,3*n) == 0
        Wh = W;
        fprintf('iter = %i, gap = %f\n', nt, mean(abs(Wh-W0)));
        if mean(abs(Wh-W0)) < 1e-3
            break;
        else
            W0 = Wh;
            rho = rho*0.95;
        end
    end
end

B = sparse(edgerow,edgecol,a2/norm(exp(Wh))*exp(Wh),p,p);
B = B+B.';
    
    
    
    