function lobjf = logobjfunc(cplmat,B,Jp)

% compute objective function
% Yu Hang, Jan, 2015, NTU

p = size(Jp,1);
n = size(cplmat,1);
[edgerow,edgecol] = find(triu(Jp,1));

l_array = zeros(n,1);

for i = 1:n
    BC = B.*sparse(edgerow,edgecol,cplmat(i,:).',p,p);
    BC = BC+BC.';
    QBC = spdiags(sum(BC,2),0,-BC);
    QBC = QBC(1:p-1,1:p-1);
    l_array(i) = logdet(QBC);
end

QB = spdiags(sum(B,2),0,-B);
QB = QB(1:p-1,1:p-1);
lobjf = mean(l_array)-logdet(QB);