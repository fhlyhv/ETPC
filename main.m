%clear
load MS8x8 Xtrain; %JPpeaks16x16_thr100 XDat;

XDat = Xtrain;
[n,~] = size(XDat);
pr = 8; % no. of rows at the bottom scale
pc = 8; % no. of columns at the bottom scale

%% marginal analysis: GEV fitting and smoothing
N = 3000;   %no. of bootstrap subsets
Id = bootstrap(XDat,N);
[L0,G0,S0,L_Var,G_Var,S_Var] = GEVPrm_bootstrap (XDat,N,Id); 
Jp = thin_membrane(pc,pr);

[Lh,alpha_u]=EM_Smth(L0,L_Var,Jp);
[Gh,alpha_g]=EM_Smth(G0,G_Var,Jp);
[Sh,alpha_s]=EM_Smth(S0,S_Var,Jp);

% compute CDF
Fx = gevcdf(XDat,repmat(Gh,n,1),repmat(Sh,n,1),repmat(Lh,n,1));

%% joint analysis using ETPC model

Fx(Fx>1-1e-8) = 1-1e-8;
Fx(Fx<1e-8) = 1e-8;

% compute the density of the pairwise copula in the ETPC model
Ks=thin_membrane(pc,pr);
[copulamatrix0, rho_mat0, nu_mat0, typ_mat0] = copuladensitymat_bic(Ks,Fx,0); % 

% compute the \Beta matrix
a1 = 50; %penalty parameter of the Lagrangian multipler
a2 = ceil(1e-1*sqrt(2*pr*(pc-1))); % |U(\Beta)|_2 = a2;
Beta0 = SGA_ssav(Ks,copulamatrix0,a1,a2);

%% missing data imputation
load MS8x8 Xtest;

Grid = ones(8);
Grid(3:6,3:6) = NaN;
XMid = find(isnan(Grid(:)));
XO = Xtest(1,:);
XO(XMid) = NaN;

[XO,lv] = SVI_ADAM(XMid,XO,Beta0,full(rho_mat0),full(nu_mat0),typ_mat0,Gh,Sh,Lh);