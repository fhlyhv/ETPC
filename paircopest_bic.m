function [i,j,rhohat,nuhat,typ,coppdf] = paircopest_bic(k,Fx,edgerow,edgecol,typi)    


i = edgerow(k);
j = edgecol(k);
U = [Fx(:,i),Fx(:,j)];

% non-Gaussian case
if typi >=1 && typi <=8
    typ = typi;
    rhohat = 0;
    nuhat = 0;
    coppdf = zeros(size(U,1),1);
    
    if typi == 1 % Gaussian
        rho_mat = copulafit('Gaussian',U);
        rhohat = rho_mat(1,2);
        coppdf = copulapdf('Gaussian',U,rhohat);
        
    elseif typi == 2 % t
        [rho_mat,nuhat] = copulafit('t',U);
        rhohat = rho_mat(1,2);
        coppdf = copulapdf('t',U,rhohat,nuhat);

    elseif typi == 3 % Gumbel
        rhohat = copulafit('Gumbel',U);
        coppdf = copulapdf('Gumbel',U,rhohat);

    elseif typi == 4 % Clayton
        rhohat = copulafit('Clayton',U);
        coppdf = copulapdf('Clayton',U,rhohat);

    elseif typi == 5 % Frank
        rhohat = copulafit('Frank',U);
        coppdf = copulapdf('Frank',U,rhohat);

    elseif typi ==6 % Galambos
        rhohat = Galamboscopulafit(U);
        coppdf = Galamboscopulapdf(U,rhohat);
        
    elseif typi == 7 % HusterReiss
        rhohat = HusterReisscopulafit(U);
        coppdf = HusterReisscopulapdf(U,rhohat);

    elseif typi == 8 % tev 
        [rhohat,nuhat] = tevcopulafit_new(U,4);
        coppdf = tevcopulapdf(U,rhohat,nuhat);
    end
    
else
    [rhohat,nuhat,typ,coppdf] = allcopulafit_bic(U,typi);
end