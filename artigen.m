clear;

%% predefine
nr = 32;
nc = 32;
p = nr*nc;
n = 1000;
K = 10000;

%% determine the underlying Gaussian process
[Ltt,Lng] = meshgrid((1:nr),(1:nc));
Ltt = Ltt(:);
Lng = Lng(:);

Phi = 8;
S = 0.01*expcor(Ltt,Lng,Phi);

%% sample max-stable process
Z = zeros(n,p);
parfor i = 1:n
    ee = exprnd(1,K,1);
    pp = cumsum(ee);
    W = mvnrnd(ones(1,p),S,K);
    Z(i,:) = max(spdiags(1./pp,0,K,K)*W);
end

%% transform the unit Frechet marginals to general GEV marginals

nc = max(Lng)+0.5;
nr = max(Ltt)+0.5;
G=0.15; %0.05+0.001*(Lng-nr/2)-0.001*(Ltt-nc/2)+0.0012*((Lng-nr/2-1).^2-1.2^2)-0.0013*(Lng-nr/2+1).*(Ltt-nc/2-1)-0.0015*((Ltt-nc/2+1).^2-1.2^2);
S=15+0.005*(Lng-nr/2)+0.01*(Ltt-nc/2)-0.008*((Lng-nr/2-0.5).^2-7.5^2)-0.005*(Lng-nr/2-2).*(Ltt-nc/2+3)+0.007*((Ltt-nc/2+0.5).^2-1.5^2);
L=30+0.02*(Lng-nr/2+2)-0.03*(Ltt-nc/2-2)+0.02*((Lng-nr/2).^2-5.5^2)-0.04*(Lng-nr/2+0.5).*(Ltt-nc/2-1.5)-0.03*((Ltt-5).^2-2.5^2);

Ga = G; %repmat(G.',n,1);
Sa = repmat(S.',n,1);
La = repmat(L.',n,1);

XDat = La + Sa./Ga.*(Z.^Ga-1);