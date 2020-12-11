function B = LowRankApp(r0,W1,S1,W2,S2,nscale,cscale1,cscale2)

% generate low rank matrix B for 2D problem
% Family -- wavelet family  p1 -- no. of rows  p2 -- no of columns  r --
% rank of each scale

% Yu Hang, NTU, Jan, 2015



B = [];

for ns = 1:nscale
    
    pc1 = cscale1(ns+1)-cscale1(ns);
    WW1 = W1(:,cscale1(ns)+1:cscale1(ns+1));
    SS1 = S1(:,cscale1(ns)+1:cscale1(ns+1));
    
    r = min(r0,pc1);
    C = spdiags(sign(binornd(1,0.5,r,pc1/r)-0.5),0:-r:1-pc1,pc1,r);
    WW1 = WW1*C;
    SS1 = SS1*C;
    
    
    pc2 = cscale2(ns+1)-cscale2(ns);
    WW2 = W2(:,cscale2(ns)+1:cscale2(ns+1));
    SS2 = S2(:,cscale2(ns)+1:cscale2(ns+1));
    r = min(r0,pc2);
    C = spdiags(sign(binornd(1,0.5,r,pc2/r)-0.5),0:-r:1-pc2,pc2,r);
    WW2 = WW2*C;
    SS2 = SS2*C;
    
    B = cat(2,B,[kron(WW2,WW1),kron(WW2,SS1),kron(SS2,WW1)]);
    
end
B = cat(2,B,kron(SS2,SS1)*sign(binornd(1,0.5,1,1)-0.5));
    


