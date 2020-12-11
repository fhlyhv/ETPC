function [valh,alphah]=EM_Smth(val0,Noise_Var,Jp)

% smoothing observed parameters val0 using EM algorithm 
% Yu Hang, Mar. 2012, NTU

p=size(Jp,1);
alpha0=0;
Rinv=spdiags(1./Noise_Var.',0,p,p);
y = sparse(val0./Noise_Var).';

for nt = 1:2000
    Kpost = alpha0*Jp+Rinv;
    x = Kpost\y;
%     L = LowRankApp(4,W1,S1,W2,S2,nscale,cscale1,cscale2);
    alphah = (p-1)/(trace(Jp/Kpost)+x'*Jp*x);  %
    if abs((alphah-alpha0)/alpha0)<1e-4 || rcond(full(alpha0*Jp+Rinv))<5e-16    %/alpha0
        break;
    else
        alpha0=alphah;
    end
    if rem(nt,100) == 0;
        nt;
    end
end

if rcond(full(alpha0*Jp+Rinv))<5e-16 || nt == 2000
    valh= sum(y)/sum(diag(Rinv))*ones(1,p);
    alphah=Inf;
else
    x=(alphah*Jp+Rinv)\Rinv*val0';
    valh=x';
end