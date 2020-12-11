% Compute the pairwise copula density matrix for every pair of connected nodes in the thin-membrane model
% Input parameters: 1) pairwise t copula matrix (testimate)
%                   2) GEV location (muGEV) 3) scale (sigmaGEV) 4) shape (xiGEV) matrix 
%                   5) realizations (x) 

%                   muGEV, sigmaGEV, and xiGEV are column vectors

% Output: 3-D pairwise copula density matrix where each 2-D matrix is
%         symmetric

% Output is matrix w as in the paper

function [coppdf, rho_mat, nu_mat, typ_mat] = copuladensitymat_bic(Jp,Fx,typi)  %

[n,p] = size(Fx);

% copulamatrix = cell(n,1); %zeros(p,p,n);

[edgerow,edgecol] = find(triu(Jp,1));
nedge = length(edgerow);
id_i = zeros(nedge,1);
id_j = zeros(nedge,1);
rhohat = zeros(nedge,1);
nuhat = zeros(nedge,1);
typ = zeros(nedge,1);
coppdf = zeros(n,nedge);

for k = 1:nedge;
    [id_i(k),id_j(k),rhohat(k),nuhat(k),typ(k),coppdf(:,k)] = paircopest_bic(k,Fx,edgerow,edgecol,typi);    
end


rho_mat = sparse(p,p);
nu_mat = sparse(p,p);
typ_mat = sparse(p,p);

rho_mat(sub2ind([p,p],id_i,id_j)) = rhohat;
rho_mat = rho_mat + rho_mat.';
nu_mat(sub2ind([p,p],id_i,id_j)) = nuhat;
nu_mat = nu_mat + nu_mat.';
typ_mat(sub2ind([p,p],id_i,id_j)) = typ;
typ_mat = typ_mat + typ_mat.';



% for i = 1:n
%     for k = 1:nedge
%         copmat = zeros(p,p);
%         copmat(id_i(k),id_j(k),:) = coppdf(i,k); 
%         copmat(id_j(k),id_i(k),:) = coppdf(i,k);
%     end
%     copulamatrix{i} = sparse(copmat);
% end



% rho_mat = sparse(rho_mat);
% tnu = sparse(tnu);

