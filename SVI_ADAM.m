function [XO,lv] = SVI_ADAM(XMid,XO,Bopt,rho_mat,nu_mat,typ_mat,Gh,Sh,Lh)

% variational inference of missing data
% Yu Hang, NTU, Jan, 2016
tic;
XMid = sort(XMid);
p = length(Gh);
pM = length(XMid);
tol = 1e-2;

% compute initial estimates as the mean of neghboring nodes
XOid = setdiff(1:p,XMid);
B = spdiags(sum(Bopt,2),0,-Bopt);
XO(XMid) = -XO(XOid)*B(XOid,XMid)/B(XMid,XMid);

Xt = XO;
XOold = XO;

% initialize SVI
m = XO(XMid).';
m0 = m;
lv = log(10*ones(pM,1));

UO = gevcdf(XO,Gh,Sh,Lh).';
ide = find(triu(Bopt,1));
nedges = length(ide);
[edgerow,edgecol] = ind2sub([p,p],ide);
Bv = Bopt(ide);
rhov = rho_mat(ide);
nuv = nu_mat(ide);
typv = typ_mat(ide);

cplv = zeros(nedges,1);
typset = [];
for i = 1:8
    if sum(typv == i) >0
       cplv(typv == i) = allcopulapdf(i,[UO(edgerow(typv==i)),UO(edgecol(typv==i))],rhov(typv==i),nuv(typv==i));
       typset = cat(2,typset,i);
    end
end

GM = Gh(XMid).';
SM = Sh(XMid).';
LM = Lh(XMid).';



BMadj = sparse(edgerow,edgecol,(1:nedges).',p,p);
XOid = setdiff(1:p,XMid);
BMadj(XOid,XOid) = 0;
ideM = BMadj(BMadj>0);
neM = length(ideM);
erM = edgerow(ideM);
ecM = edgecol(ideM);
cplvm = zeros(neM,1);
E = sparse(neM,p);
for i = 1:neM
    if max(erM(i),ecM(i)) < p
        E(i,erM(i)) = 1;
        E(i,ecM(i)) = -1;
    else
        E(i,min(erM(i),ecM(i))) = 1;
    end
end
BeM = Bv(ideM);
typeM = typv(ideM);
rhoeM = rhov(ideM);
nueM = nuv(ideM);

BMadj = sparse(erM,ecM,(1:neM).',p,p);
BMadj = BMadj+BMadj.';
BMadj(setdiff(1:p,XMid),:) = 0;
[vcM,vrM] = find(BMadj.'>0);  % must be indexed by rows
nvM = length(vrM);
pcpu = zeros(nvM,1);
S = sparse(pM,nvM);  % selection matrix
k = 0;
for i = 1:pM
    ine = sum(sign(BMadj(XMid(i),:)));
    S(i,k+1:k+ine) = 1;
    k = k+ine;
end
idvM = BMadj(sub2ind([p,p],vrM,vcM));
typvM = typeM(idvM);
rhovM = rhoeM(idvM);
nuvM = nueM(idvM);

% determine whetehr to use low-rank approximation or not
K = Bopt + diag(sum(Bopt)) + speye(p);
c = graphcolor_irregular(K^4);
if neM > max(c)
    flag1 = 1;
else
    flag1 = 0;
end

% generate quasi random numbers
if pM>6
    urnd = sobolset(pM,'Skip',1e3,'Leap',1e2);
    urnd = scramble(urnd,'MatousekAffineOwen');
    urnd = qrandstream(urnd);
else
    urnd = haltonset(pM,'Skip',1e3,'Leap',1e2);
    urnd = scramble(urnd,'RR2');
    urnd = qrandstream(urnd);
end



gx = 0;
pppv = 0;
dx1 = 0;
dv1 = 0;
dx2 = 0;
dv2 = 0;
tau1 = 0.9;
tau2 = 0.999;
dtau = sqrt(1-tau2)/(1-tau1);
eta = 0.01;   % Pls reduce the step size if the algorithm diverges
v0 = 0;
eps = 1e-6*(1-tau2);


for nt = 1:1e7
    v = exp(lv);
    z = norminv(qrand(urnd,1)).';
    Xtmp = v.*z + m;    
    a0 = 1+GM./SM.*(Xtmp-LM);
    if sum(a0(:)<=0) > 0
        continue;
    end
    
    
    
    utmp = gevcdf(Xtmp,GM,SM,LM);
    utmp(utmp==0) = 1e-6;
    UO(XMid) = utmp;
    for i = typset
        if any(typeM==i)
            cplvm(typeM==i) = allcopulapdf(i,[UO(erM(typeM==i)),UO(ecM(typeM==i))],rhoeM(typeM==i),nueM(typeM==i));
        end
        if any(typvM==i)
            Uu = min(UO(vrM(typvM==i))+1e-4,1-1e-6);
            Ul = max(UO(vrM(typvM==i))-1e-4,1e-6);
            pcpu(typvM==i) = (allcopulapdf(i,[Uu,UO(vcM(typvM==i))],rhovM(typvM==i),nuvM(typvM==i))-...
                allcopulapdf(i,[Ul,UO(vcM(typvM==i))],rhovM(typvM==i),nuvM(typvM==i)))./(Uu-Ul);
        end
    end

    cplv(ideM) = cplvm;
    BC = sparse(edgerow,edgecol,Bv.*cplv,p,p);
    BC = BC+BC.';
    QBC = spdiags(sum(BC,2),0,-BC);
    QBC = [QBC(1:p-1,1:p-1),sparse(p-1,1);sparse(1,p-1),1];
    %[~,q] = chol(QBC);
    q = condest(QBC);
    if q > 1e6 || isnan(q)
        continue;
    end
    
    if nt > 2.6e4
        q;
    end
    
    if flag1 == 1
        L = sparse((1:p).',c,sign(randi([0,1],p,1)-0.5)); 
        EL = E*L;
        RC = QBC\L;
        pppc = BeM.*sum((E*RC).*EL,2);
    else
        pppc = BeM.*sum((E/QBC).*E,2);
    end
    pppu = S*(pcpu.*pppc(idvM));
    
    pfpx = -(1+GM)./SM./a0 + a0.^(-1./GM-1)./SM;
    pupx = gevpdf(Xtmp,GM,SM,LM);
    
    pppx = pppu.*pupx+pfpx;
    pppv = pppx.*z.*v+1;
    

    dx1 = tau1*dx1+(1-tau1)*pppx;
    dx2 = tau2*dx2+(1-tau2)*pppx.^2;
    
    dv1 = tau1*dv1+(1-tau1)*pppv;
    dv2 = tau2*dv2+(1-tau2)*pppv.^2;

    m = m + eta*dtau.*dx1./sqrt(dx2+eps);
    lv = lv + eta*dtau.*dv1./sqrt(dv2+eps);

    
    if rem(nt,1e4) == 0
        nt;
    end
    
    if rem(nt,1e3) == 0
        mh = m;
        fprintf('iter = %i, gap = %f\n', nt, max(abs(mh-m0)));
        if   max(abs(mh-m0)) < tol
            t = toc
            break;
        else
            m0 = mh;
            eta = eta*0.99;
        end
        
    end
end
XO(XMid) = m;            