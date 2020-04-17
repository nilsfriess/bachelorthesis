clf

N = 200;
a = sprandsym(N,0.3);

% a = load("fe_matrix.mat");
% a = a.Ahat;
% [N,~] = size(a);

[V,D] = eigs(a,N);
D = diag(D); % D is diagonal matrix, extract eigenvalues and store back

% Create random linear combination of eigenvectors of A
% where one of the interior eigenvectors is weighted heavier
% than the rest (we later expect convergence to that very evec)
targetIndex = randi(N);
%weights = rand(N,1);
%weights(targetIndex) = 10;
targetV = V(:,targetIndex);
targetE = D(targetIndex);

disp(['Target eigenvalue: ', num2str(targetE)]);

angle = 45;
v = random_in_cone(targetV, angle);

disp(['Residual = ', num2str(norm(a*v - rayleighquotient(a,v)*v))]);

disp(['Angle between v and target (deg): ', num2str(rad2deg(acos(v'* targetV)))]);
disp(['Angle between v and target (rad): ', num2str((acos(v'* targetV)))]);
disp(' ');

% sqrt(1 - (v'*targetV)^2)

gamma = 1;
atilde = a - gamma*1i*(speye(N) - v*v');
[Vtilde, Dtilde] = eigs(atilde, N);
Dtilde = diag(Dtilde); % extract eigenvalues

% Test classic RQI
tic; [x1, e1, its1, e_its1, v_its1, res1] = classic_rqi(a, v, inf); toc;
disp('Result of classic RQI');
disp(['Computed eigenvalue: ', num2str(e1), ' (', num2str(its1), ' Iterations)']);
disp(' ');

q = 1;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)

% Test complex RQI with gamma = res
tic; [x2, e2, its2, e_its2, v_its2, res2] = complex_rqi2(a, v, inf, inf); toc;
disp('Result of complex RQI (gamma = res)');
disp(['Computed eigenvalue: ', num2str(e2), ' (', num2str(its2), ' Iterations)']);
if abs(targetE - e2) < 10e-9
    disp("Complex RQI succeded");
else
    disp("Complex RQI did not succeed");
end

q = 1;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)
q = 2;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)
q = 3;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)

q = 1;
diff2 = (real(e_its2(2:end)) - targetE) ./ (real(e_its2(1:end-1)) - targetE).^q 

disp(' ');

% Test complex RQI with gamma = res^2
tic; [x3, e3, its3, e_its3, v_its3, res3] = complex_rqi2square(a, v, inf, inf); toc;
disp('Result of complex RQI (gamma = res^2)');
disp(['Computed eigenvalue: ', num2str(e3), ' (', num2str(its3), ' Iterations)']);
if abs(targetE - e3) < 10e-9
    disp("Complex RQI (square) succeded");
else
    disp("Complex RQI (square) did not succeed");
end

q = 1;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)
q = 2;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)
q = 3;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)

disp(' ');


% Test complex RQI with gamma adaptively
tic; [x4, e4, its4, e_its4, v_its4, res4] = complex_rqi2combined(a, v, inf, inf); toc;
disp('Result of complex RQI (Combined)');
disp(['Computed eigenvalue: ', num2str(e4), ' (', num2str(its4), ' Iterations)']);
if abs(targetE - e4) < 10e-9
    disp("Complex RQI (combined) succeded");
else
    disp("Complex RQI (combined) did not succeed");
end

q = 1;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)
q = 2;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)
q = 3;
diff = res1(2:end) ./ (res1(1:end-1) .^ q)

% Plot residuals
semilogy(1:its1+1, res1, '-o'); hold on
semilogy(1:its2+1, res2, '-x')
semilogy(1:its3+1, res3, '-s')
semilogy(1:its4+1, res4, '-^')
pbaspect([1 1 1])

xlabel("Iteration"); ylabel("Residual norm");

axis([0 (max([its1, its2, its3, its4]) + 2) 10^(-17) 10^2])

legend({'Classic RQI', ...
        'Complex RQI (\gamma^{(k)} = r)', ...
        'Complex RQI (\gamma^{(k)} = r^2)', ...
        'Complex RQI (\gamma^{(k)} adaptive)'}, ...
        'Location', 'southwest');
%export_fig 'crqi_residuals.eps' -transparent