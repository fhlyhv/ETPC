load resultarti32x32_11 Xtest Bopt0 rho_mat0 nu_mat0 typ_mat0 ... %Bopt1 rho_mat1 nu_mat1 typ_mat1 ...
    Gh Sh Lh W1 S1 W2 S2 nscale cscale1 cscale2 s2 SS

nt = 1;  % no. of trials
[nq,p] = size(Xtest);

MAE0_arrayr = zeros(nt,nq);
MAE1_arrayr = zeros(nt,nq);
MAEg_arrayr = zeros(nt,nq);

Xest0_arrayr = zeros(nt,p,nq);
Xest1_arrayr = zeros(nt,p,nq);
Xestg_arrayr = zeros(nt,p,nq);

k =40;
for i = 1:100
    idn = randperm(p);
    idn = sort(idn(1:k));
    for j = 3:nq
        Xest0_arrayr(i,:,j) = SVI(idn,Xtest(j,:),Bopt0,rho_mat0,nu_mat0,typ_mat0,Gh,Sh,Lh,W1,S1,W2,S2,nscale,cscale1,cscale2,s2); %ICM(XMid{i},Xtest(j,:),Gh,Sh,Lh,Bopt0,typ_mat0,rho_mat0,nu_mat0,s2,1e-2);
        Xest1_arrayr(i,:,j) = SVI(idn,Xtest(j,:),Bopt1,rho_mat1,nu_mat1,typ_mat1,Gh,Sh,Lh,W1,S1,W2,S2,nscale,cscale1,cscale2,s2);%ICM(XMid{i},Xtest(j,:),Gh,Sh,Lh,Bopt1,typ_mat1,rho_mat1,nu_mat1,s2,1e-2);
        Xestg_arrayr(i,:,j) = CGGM_impt(Xtest(j,:),idn,Gh,Sh,Lh,inv(cov_normalize(SS)));
        MAE0_arrayr(i,j) = sum(abs(Xest0_arrayr(i,:,j)-Xtest(j,:)))/k; %i^2;
        MAE1_arrayr(i,j) = sum(abs(Xest1_arrayr(i,:,j)-Xtest(j,:)))/k; %i^2;
        MAEg_arrayr(i,j) = sum(abs(Xestg_arrayr(i,:,j)-Xtest(j,:)))/k; %i^2;
    end
%     save pr;
end