function Id = bootstrap(XDat,N)

n=size(XDat,1);
Id=zeros(n,N);

for i=1:N 
    I=floor(rand(n,1)*n)+1;
    I(I==n+1)=n;
    Id(:,i)=I;
end


% [n,p]=size(XDat);
% l=floor(n/2);
% %X=zeros(n,p,N);
% Id=zeros(l,N);
% 
% for i=1:2:N-1
%     I=randperm(n)';
%     Id(:,i)=I(1:l);
%     Id(:,i+1)=I(l+1:2*l);
% end