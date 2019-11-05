clear;

topright = 1.1;
[nodes,elements] = initial_mesh(topright);

for i = 1:2
    [nodes,elements] = refine(nodes,elements);
end

% First assemble the element stiffness matrices
elt_matrices = elt_stiffness(elements, nodes);

% Next assemble the global stiffness matrix 
% (including all boundary nodes)
Ahat = global_stiffness(elt_matrices, elements, nodes);

[v,d] = eigs(Ahat, 20);
eigvals = diag(d);
eval10 = eigvals(10);
eval11 = eigvals(11);

disp(['10th eigenvalue = ', num2str(eval10)])
disp(['11th eigenvalue = ', num2str(eval11)]) 

sigma = eval11 + 0.1;

disp(' ')
disp('TEST #1')
disp(['Using classic RQI with shift s = ', num2str(sigma)])

evec10 = v(:,10);
evec11 = v(:,11);

v0 = evec10;
v0(10) = v0(10) + 0.2;
v0 = v0 / norm(v0);

rad2deg(acos(evec10'*v0))
rad2deg(acos(evec11'*v0))