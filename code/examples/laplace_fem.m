% Now, create model for slightly perturbed problem
modelP = createpde(1);
% Geometry for "perturbed" unit square
upperright = 2;
gd = [3;4;0;1;1;0;1;upperright;0;0];
ns = 'S';
sf = 'S';
g = decsg(gd,ns,sf');
geometryFromEdges(modelP,g);
% Laplace eigenvalues
specifyCoefficients(modelP,'m',0,'d',1,'c',1,'a',0,'f',1);

% Apply Neumann and Dirichlet BC's
applyBoundaryCondition(modelP, 'neumann', 'Edge', 1, 'g', 1);
applyBoundaryCondition(modelP, 'edge', [2,3,4], 'u', 0);

generateMesh(modelP,'Hmax',0.05);


X = modelP.Mesh.Nodes(1,:);
Y = modelP.Mesh.Nodes(2,:);

k = 2; l = 2;
laplace_eigval_exact(k,l);
u = laplace_eigfun_exact(k,l,X,Y);
u = u / norm(u);

%figure(1);
%pdeplot(modelP, 'XYData', u, 'ZData', u);

FEM = assembleFEMatrices(modelP, 'nullspace');
K = FEM.Kc; B = FEM.B; M = FEM.M;
[V,E] = eigs(K,M,1,59);
E = diag(E);
V = B*V;

vEigs = V(:,1);
E

result = solvepdeeig(modelP,[55,65]);
result.Eigenvalues

v = result.Eigenvectors(:,1);
v = v / norm(v);

norm(v - vEigs)

%figure(2);
%pdeplot(modelP, 'XYData', v, 'ZData', v);

size(K)
utest = B'*u';
utest = utest / norm(utest);
size(utest)

vtest = B'*v;
vtest = vtest / norm(vtest);

angle = rad2deg(acos(abs(vtest'*utest)))

% Test of CRQI
[x, e, its, e_its, x_its, res] = complex_rqi2combined_general(K, M, utest, inf, inf, 10e-11);

e



