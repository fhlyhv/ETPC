function c = graphcolor_regular(pr,pc,hcl)

% graph coloring for regular grid using 18 colors
% Yu Hang, NTU, Jan, 2015
% pr -- no. of rows; pc -- no. of columns; hcl -- half correlation length


c = zeros(pr,pc);
cl = 2*hcl;

for i = 1:hcl
    for nc = 1:cl
        c(i:cl:pr,nc:cl:pc) = (i-1)*cl+nc;
        if nc <= hcl
            c(i+hcl:cl:pr,nc+hcl:cl:pc) = (i-1)*cl+nc;
        else
            c(i+hcl:cl:pr,nc-hcl:cl:pc) = (i-1)*cl+nc;
        end
    end
end

c = c.';
c = c(:);