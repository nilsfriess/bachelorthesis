clf

test = 4;
size = 100;
weight_target = 10;
gamma = 1;

[a,N] = test_matrix(test,size);

[V,D] = eigs(a,N);
D = diag(D); % D is diagonal matrix, extract eigenvalues and store back

% Create random linear combination of eigenvectors of A
% where one of the interior eigenvectors is weighted heavier
% than the rest (we later expect convergence to that very evec)
targetIndex = 45;
weights = 8*rand(N,1);
weights(targetIndex) = weights(targetIndex) + weight_target;
targetV = V(:,targetIndex);
targetE = D(targetIndex);

v = randn(N,1);
v = v + weight_target * targetV;
v = v / norm(v);
fprintf("Matrix size: %d\n\n", N);

disp(['Target eigenvalue: ', num2str(targetE)]);


disp(['Residual = ', num2str(norm(a*v - rayleighquotient(a,v)*v))]);

disp(['Angle between v and target (deg): ', num2str(rad2deg(acos(v'* targetV)))]);
disp(['Angle between v and target (rad): ', num2str((acos(v'* targetV)))]);
disp(' ');

Atilde = a - gamma*1i*(speye(N) - v*v');

Dtilde = eigs(Atilde, N);

plot(complex(D), 'or', 'Markersize', 5, 'MarkerFaceColor', 'r');
hold on;

plot(Dtilde, 'ob', 'Markersize', 5, 'MarkerFaceColor', 'b');

pbaspect([3,1,1]);
axis([0,16,-1.3,0.3]);

legend(["Spectrum of $A$", "Spectrum of $\tilde{A}$"], 'Interpreter', 'latex')

xlabel("Real");
ylabel("Imaginary");