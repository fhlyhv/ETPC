function [L0,G0,S0,L_Var,G_Var,S_Var] = GEVPrm_bootstrap (XDat,N,Id)

%compute ML estimate and bootstrap variance of each parameter
%Yu Hang, Mar. 2012, NTU

%% predifine
[n,p]=size(XDat);
L_bt=zeros(N,p);  
G_bt=zeros(N,p);
S_bt=zeros(N,p);


%% GEV fit
% L0=zeros(1,p);  %location
% G0=zeros(1,p);  %shape
% S0=zeros(1,p);  %scale


% parfor i=1:p
%     [G0(i),S0(i),L0(i)]=gevfit_med_mex(XDat(:,i));
% end

[L0,G0,S0]=X_gevfit(XDat);

%% Bootstrap variacne


parfor i=1:N
%     X = gevrnd(repmat(G0,n,1),repmat(S0,n,1),repmat(L0,n,1));
    [L_bt(i,:),G_bt(i,:),S_bt(i,:)]=X_gevfit(XDat(Id(:,i),:));  %  
end


L_Var=var(L_bt);
G_Var=var(G_bt);
S_Var=var(S_bt);