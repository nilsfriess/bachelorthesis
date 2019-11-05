function [newnodes,newelements] = refine(nodes,elements)

[nnodes,m]  = size(nodes);
[nelts,m]   = size(elements);

% cycle over the elements; find coordinates of new nodes
% and add them to the array nodes

newnodes = nodes;
inode = 1;
edgetable = zeros(0,2);  % Create an empty matrix with two columns

for ie = 1:nelts

  % First of all look at all possible edges of element ie
  % if they have been visited before, a new node has already 
  % been created at the midpoint of the edge. The function
  % "visited" returns the (global) number of this node
    
  for edge = 1:3    

    n1 = edge;
    n2 = mod(edge,3)+1;
    
    newloc(edge) = visited(edgetable,...
                           elements(ie,n1),elements(ie,n2));

    % If the edge has not been visited before create a new
    % node at the midpoint of the edge and mark it as 
    % Dirichlet to begin with. If it is visited again then 
    % it must be an interior edge.
    
    if (newloc(edge) == 0)
        
        newloc(edge) = inode;
        inode = inode+1;
        
        edgetable = [edgetable;elements(ie,n1),elements(ie,n2)];
                 
        newnodes = [newnodes;(nodes(elements(ie,n1),1:2) +...
                    nodes(elements(ie,n2),1:2))/2,1];
                
    else
        newnodes(nnodes+newloc(edge),3) = 0;
    end
    
  end 
  
  newloc = newloc + nnodes;
        
  % Once all the edges of element ie have been refined create 
  % the four refined elements (care has to be taken to make 
  % sure the local numbering is anticlockwise!)
  
  for tau = 1:3
      
      newelements(4*(ie-1)+tau,1:3) = [elements(ie,tau),...
          newloc(tau),newloc(mod(tau+1,3)+1)];
 
  end   
  
  newelements(4*ie,1:3) = newloc(1:3);

end
