% Now, create model for slightly perturbed problem
modelP = createpde(1);
% Geometry for "perturbed" unit square
gd = [3;4; 0; 0.9;...
           1; 0;...
           1; 1.1;...
           0; 0.1];
ns = 'S';
sf = 'S';
g = decsg(gd,ns,sf');
geometryFromEdges(modelP,g);

% Setup Laplace eigenvalue problem
specifyCoefficients(modelP,'m',0,'d',1,'c',1,'a',0,'f',1);

% Apply Neumann and Dirichlet BC's (three fixed and one free side)
applyBoundaryCondition(modelP, 'neumann', 'Edge', 1, 'g', 1);
applyBoundaryCondition(modelP, 'edge', [2,3,4], 'u', 0);

generateMesh(modelP,'Hmax',0.04);

X = modelP.Mesh.Nodes(1,:);
Y = modelP.Mesh.Nodes(2,:);

k = 2; l = 2;
e = laplace_eigval_exact(k,l);
u = laplace_eigfun_exact(k,l,X,Y);

% for i = 1:length(u)
%     if u(i) > 0
%         u(i) = 1;
%     else
%         u(i) = -1;
%     end
% end

u = u' / norm(u);

subplot(3,1,1);
pdeplot(modelP, 'XYData', u, 'ZData', u);
title("Exact eigenfunction of Laplace on unit square")

FEM = assembleFEMatrices(modelP, 'nullspace');
K = FEM.Kc; B = FEM.B; M = FEM.M;
[V,E] = eigs(K,M,1,e);
E = diag(E);
V = B*V;
size(K)

vEigs = V(:,1);
 
result = solvepdeeig(modelP,[e-10,e+10]);
fprintf("Target eigenvalue (pdeeig): %.5e\n", result.Eigenvalues(1));
fprintf("Target eigenvalue (eigs):   %.5e\n", E(1));

v = result.Eigenvectors(:,1);
v = v / norm(v);

subplot(3,1,2);
pdeplot(modelP, 'XYData', v, 'ZData', v);
title("Target eigenfunction")

if (norm(v - vEigs) - 2 < 10e-8)
   v = -v; 
end

fprintf("Difference between pdeeig result and eig result: %.5e\n", norm(v - vEigs));

fprintf("Angle between exact and guess: %.2f\n", rad2deg(acos(abs(u'*v))));

% Remove Boundary conditions from vectors
utest = B'*u;
utest = utest / norm(utest);

vtest = B'*v;
vtest = vtest / norm(vtest);

% Test of CRQI
[x, e, its, e_its, x_its, res] = complex_rqi2combined_general(K, M, utest, inf, inf, 10e-11);

x = B*x;

% Right eigenfunction but wrong sign
norm(x - v)
if (norm(x - v) - 2) < 10e-8
    x = -x;
end

subplot(3,1,3);
pdeplot(modelP, 'XYData', x, 'ZData', x);
title("Computed eigenfunction")




