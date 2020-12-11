function Strue=expcor(Ltt,Lng,Phi)

p=length(Ltt);
Strue=exp(-((repmat(Ltt,1,p)-repmat(Ltt.',p,1)).^2+(repmat(Lng,1,p)-repmat(Lng.',p,1)).^2)/Phi);



% 
% for i=1:p-1
%     for j=i+1:p
%         Strue(i,j)=exp(-sum(([Ltt(i),Lng(i)]-[Ltt(j),Lng(j)]).^2)/Phi);  %sqrt
%         Strue(j,i)=Strue(i,j);
%     end
% end