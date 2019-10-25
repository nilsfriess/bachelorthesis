function Ahat = global_stiffness(elt_matrices,elements,nodes)

% function which assembles global stiffness matrix

[nelts,m] = size(elements);
[nnodes,m] = size(nodes);

% Create an empty sparse matrix with nnodes columns and rows

Ahat = sparse(nnodes,nnodes);

% Assemble the global stiffness matrix from the element matrices
% including all the boundary nodes

for ie=1:nelts

    iglob = elements(ie,:);
    Ahat(iglob,iglob) = Ahat(iglob,iglob) + elt_matrices{ie};
   
end
