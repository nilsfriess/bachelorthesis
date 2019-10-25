function f = rhs(func,N0,elements,nodes);

% assembles the load vector arising from the function specified 
% in 'func' on the mesh specified by elements and nodes
% There is an entry of the output vector f for each node
% in the list N0

[nelts,m]  = size(elements);  % Finds how many elements there are
[nnodes,m] = size(nodes);     % Finds how many nodes there are

fhat = zeros(nnodes,1);       % initialise with a vector of zeros

for ie = 1:nelts
    
    X = nodes(elements(ie,:)',1:2);   % Coordinates of 3 nodes of 
                                      % element ie

    area = mu(X);                     % Area of element ie
        
    for p = 1:3
        
        i = elements(ie,p);
        fhat(i) = fhat(i) + feval(func,X(p,1),X(p,2)) * area/3;
            
    end
   
end

f = fhat(N0);
