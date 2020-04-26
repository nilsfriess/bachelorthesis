a = sprandsym(10,0.5);
%a =  diag([4,9,22,11,44,1,5,2,8,12]);

weights = randn(10,1);
weights(5) = 10;

[V,D] = eigs(a,10,'sr');
D = diag(D);
v5 = V(:,5);

v = V*weights;
v = v/ norm(v);

atilde = a - 1i*(eye(10) - v*v');
[Vtilde,Dtilde] = eigs(atilde,10,'sr');
Dtilde = diag(Dtilde);

atilde0 = a - 1i*(eye(10) - v5*v5');
atilde1 = -1i*(v5*v5' - v*v');

[Vtilde0,Dtilde0] = eigs(atilde0,10,'sr');
Dtilde0 = diag(Dtilde0);

[Vtilde1,Dtilde1] = eigs(atilde1,10,'sr');
Dtilde1 = diag(Dtilde1);
