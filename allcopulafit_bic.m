function [rho,nu,typ,PDF] = allcopulafit_bic(U,typi)

% typ: 
%      elliptical
%      1 -- Gaussian;
%      2 -- t;
%      Archimedemean
%      3 -- Gumbel;
%      4 -- Clayton;
%      5 -- Frank;
%      Extreme value
%      6 -- Galambos;
%      7 -- HusterReiss;
%      8 -- tev

% typi
%      0 -- min BIC
%      9 -- tail dependence
%      10-- max stable
U(U>=1) = 1-1e-16;
U(U<=0) = 1e-16;
n = length(U);
bic = zeros(1,8);
rho_array = zeros(1,8);
nu_array = zeros(1,8);
PDF_array = zeros(n,8);

% Gaussian
rho_mat = multicopulafit('Gaussian',U);
rho_array(1) = rho_mat(1,2);
if rho_array(1)>=1
    rho_array(1) = 1-1e-6;
elseif rho_array(1)<=-1
    rho_array(1) = -1+1e-6;
end

PDF_array(:,1) = copulapdf('Gaussian',U,rho_array(1));
bic(1) = -2*sum(log(PDF_array(:,1))) + log(n);

% t
[rho_mat,nu_array(2)] = multicopulafit('t',U);
rho_array(2) = rho_mat(1,2);

if rho_array(2)>=1
    rho_array(2) = 1-1e-6;
elseif rho_array(2)<=-1
    rho_array(2) = -1+1e-6;
end
PDF_array(:,2) = copulapdf('t',U,rho_array(2),nu_array(2));
bic(2) = -2*sum(log(PDF_array(:,2))) + 2*log(n);

% Gumbel
rho_array(3) = multicopulafit('Gumbel',U);
if rho_array(3) <=1
    PDF_array(:,3) = ones(n,1);
else
    PDF_array(:,3) = copulapdf('Gumbel',U,rho_array(3));
end
bic(3) = -2*sum(log(PDF_array(:,3))) + log(n);

% Clayton
rho_array(4) = multicopulafit('Clayton',U);
if rho_array(4) <=0
    PDF_array(:,4) = ones(n,1);
else
    PDF_array(:,4) = copulapdf('Clayton',U,rho_array(4));
end
bic(4) = -2*sum(log(PDF_array(:,4))) + log(n);

% Frank
rho_array(5) = multicopulafit('Frank',U);
if rho_array(5) <=0
    PDF_array(:,5) = ones(n,1);
else
    PDF_array(:,5) = copulapdf('Frank',U,rho_array(5));
end
bic(5) = -2*sum(log(PDF_array(:,5))) + log(n);

% Galambos
rho_array(6) = Galamboscopulafit(U);
PDF_array(:,6) = Galamboscopulapdf(U,rho_array(6));
bic(6) = -2*sum(log(PDF_array(:,6))) + log(n);

% HusterReiss
rho_array(7) = HusterReisscopulafit(U);
PDF_array(:,7) = HusterReisscopulapdf(U,rho_array(7));
bic(7) = -2*sum(log(PDF_array(:,7))) + log(n);

% tev 
[rho_array(8),nu_array(8)] = tevcopulafit_new(U,4);
PDF_array(:,8) = tevcopulapdf(U,rho_array(8),nu_array(8));
bic(8) = -2*sum(log(PDF_array(:,8))) + 2*log(n);


if typi ==0
    [~,typ]= min(bic);
    rho = rho_array(typ);
    nu = nu_array(typ);
    PDF = PDF_array(:,typ);
            
elseif typi >=1 && typi <=8
    typ = typi;
    rho = rho_array(typ);
    nu = nu_array(typ);
    PDF = PDF_array(:,typ);
elseif typi ==9
    typ_array = [2,3,6,7,8];
    [~,q] = min(bic(typ_array));
    typ = typ_array(q);
    rho = rho_array(typ);
    nu = nu_array(typ);
    PDF = PDF_array(:,typ);
elseif typi ==10
    typ_array = [3,6,7,8];
    [~,q] = min(bic(typ_array));
    typ = typ_array(q);
    rho = rho_array(typ);
    nu = nu_array(typ);
    PDF = PDF_array(:,typ);
end

