function [rho,nu,typ,PDF] = allcopulafit(U,typi)

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
%      0 -- min KL div
%      9 -- tail dependence
%      10-- max stable

n = length(U);
kldiv = zeros(1,8);
rho_array = zeros(1,8);
nu_array = zeros(1,8);
PDF_array = zeros(n,8);

% Gaussian
rho_mat = copulafit('Gaussian',U);
rho_array(1) = rho_mat(1,2);
PDF_array(:,1) = copulapdf('Gaussian',U,rho_array(1));
kldiv(1) = -sum(log(PDF_array(:,1)))/n;

% t
[rho_mat,nu_array(2)] = copulafit('t',U);
rho_array(2) = rho_mat(1,2);
PDF_array(:,2) = copulapdf('t',U,rho_array(2),nu_array(2));
kldiv(2) = -sum(log(PDF_array(:,2)))/n;

% Gumbel
rho_array(3) = copulafit('Gumbel',U);
PDF_array(:,3) = copulapdf('Gumbel',U,rho_array(3));
kldiv(3) = -sum(log(PDF_array(:,3)))/n;

% Clayton
rho_array(4) = copulafit('Clayton',U);
PDF_array(:,4) = copulapdf('Clayton',U,rho_array(4));
kldiv(4) = -sum(log(PDF_array(:,4)))/n;

% Frank
rho_array(5) = copulafit('Frank',U);
PDF_array(:,5) = copulapdf('Frank',U,rho_array(5));
kldiv(5) = -sum(log(PDF_array(:,5)))/n;

% Galambos
rho_array(6) = Galamboscopulafit(U);
PDF_array(:,6) = Galamboscopulapdf(U,rho_array(6));
kldiv(6) = -sum(log(PDF_array(:,6)))/n;

% HusterReiss
rho_array(7) = HusterReisscopulafit(U);
PDF_array(:,7) = HusterReisscopulapdf(U,rho_array(7));
kldiv(7) = -sum(log(PDF_array(:,7)))/n;

% tev 
[rho_array(8),nu_array(8)] = tevcopulafit_new(U,4);
PDF_array(:,8) = tevcopulapdf(U,rho_array(8),nu_array(8));
kldiv(8) = -sum(log(PDF_array(:,8)))/n;


if typi ==0
    [~,typ]= min(kldiv);
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
    [~,q] = min(kldiv(typ_array));
    typ = typ_array(q);
    rho = rho_array(typ);
    nu = nu_array(typ);
    PDF = PDF_array(:,typ);
elseif typi ==10
    typ_array = [3,6,7,8];
    [~,q] = min(kldiv(typ_array));
    typ = typ_array(q);
    rho = rho_array(typ);
    nu = nu_array(typ);
    PDF = PDF_array(:,typ);
end

