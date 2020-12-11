function c = graphcolor_irregular(K)

% graph coloring using greedy multicoloring algorithm
% Yu Hang, NTU, Mar, 2013


p = size(K,1);
c = zeros(p,1);

for i = 1:p
    if i == 1
        c(i) = 1;
    else
        cold = c(K(1:i,i)~=0);
        c(i) = min(setdiff(1:max(cold)+1,cold));
    end
end