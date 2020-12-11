function [L,G,S]=X_gevfit(XDat)


%% predifine 
p=size(XDat,2);
L=zeros(1,p);  %location
G=zeros(1,p);  %shape
S=zeros(1,p);  %scale
% OPTIONS = statset('MaxIter',1000,'MaxFunEvals',1000);

%% GEV fit
for i=1:p
    [G(i),S(i),L(i)]=gevfit_pwm(XDat(:,i));
%     Prm=gevfit(XDat(:,i),[],OPTIONS);
%     G(i)=Prm(1);
%     S(i)=Prm(2);
%     L(i)=Prm(3);
end