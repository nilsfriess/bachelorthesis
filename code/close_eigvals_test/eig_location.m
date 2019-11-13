topright = 1.02;
[nodes,elements] = initial_mesh(topright);

for i = 1:1
    [nodes,elements] = refine(nodes,elements);
end

% First assemble the element stiffness matrices
elt_matrices = elt_stiffness(elements, nodes);

% Next assemble the global stiffness matrix 
% (including all boundary nodes)
a = global_stiffness(elt_matrices, elements, nodes);

[r,~] = size(a);
a = a + eye(r);

disp("Eigenvalues of A");
[V,D] = eigs(a, r);
diag(D)

gamma = 2;
lambda = 0;

v = V(:,5) * 2 + V(:,1) + V(:,3) + V(:,9) + V(:,2);
v = v / norm(v);

sig = lambda + 1i*gamma;
b = a - sig * eye(r);

disp("Eigenvalues of A shifted");
[~,D] = eigs(b, r);
diag(D)

%lambda = 0;

disp("Eigenvalues of A with orthognoal projection");
c = (a-lambda*eye(r)) - gamma*1i*(eye(r)-v*v');
[~,D] = eigs(c, r);
diag(D)

w = V(:,5)

w'*b^-1*w
w'*c^-1*w
