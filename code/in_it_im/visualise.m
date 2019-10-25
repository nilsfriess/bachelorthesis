function visualise(nodes,elements)

[nelts,m] = size(elements);   % Finds the number of elements

% visualises the mesh specified by the data in nodes and elements

hold off 

for ie = 1:nelts

    plot(nodes(elements(ie,1:2)',1),nodes(elements(ie,1:2)',2),...
         'k-','LineWidth',2);
    hold on
    plot(nodes(elements(ie,2:3)',1),nodes(elements(ie,2:3)',2),...
         'k-','LineWidth',2);
    plot(nodes(elements(ie,1:2:3)',1),nodes(elements(ie,1:2:3)',2),...
         'k-','LineWidth',2);

end 

% Now plot the nodes on top: use a magenta dot for Dirichlet nodes 
% and a blue dot for non-Dirichlet nodes

% Find Dirichlet nodes 
ND = find(nodes(:,3)==1);
% Plot them 
plot(nodes(ND,1),nodes(ND,2),'m.','Markersize',10)

% Find non-Dirichlet nodes
N0 = find(nodes(:,3)~=1);
% Plot them 
plot(nodes(N0,1),nodes(N0,2),'b.','Markersize',10)

hold off
