




function [elt_matrices] = elt_stiffness(elements,nodes)

% This function computes the element stiffness matrices given 
% a list of nodes and an element table. The output is a cell 
% array with an entry for each element of the mesh. Each "cell" 
% is the corresponding 3x3 element stiffness matrix.

[nelts,m] = size(elements);  % Finds how many elements there are

elt_matrices = cell(1,nelts);  

for ie = 1:nelts

    elt_matrices{ie} = zeros(3,3);

    X = nodes(elements(ie,:)',1:2);   % Coordinates of 3 nodes of 
                                      % element ie

    area = mu(X);                % Area of element ie using the 
                                 % determinant formula from lectures

    E = [X(2,:)-X(3,:); X(3,:)-X(1,:); X(1,:)-X(2,:)];  
                                 % difference vectors

    % Now we implement the formula given in lectures for the 
    % element stiffness matrices corresponding to the Laplacian. 
    
    Atau = [[E(1,:)*E(1,:)',E(1,:)*E(2,:)',E(1,:)*E(3,:)'];...
            [E(1,:)*E(2,:)',E(2,:)*E(2,:)',E(2,:)*E(3,:)'];...
            [E(1,:)*E(3,:)',E(2,:)*E(3,:)',E(3,:)*E(3,:)']];

    elt_matrices{ie} = Atau/(4*area);

end