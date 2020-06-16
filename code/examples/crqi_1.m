clf

test = 4;
size = 1000;
weight_target = 20;

[a,N] = test_matrix(test,size);

[V,D] = eigs(a,N);
D = diag(D); % D is diagonal matrix, extract eigenvalues and store back

% Create random linear combination of eigenvectors of A
% where one of the interior eigenvectors is weighted heavier
% than the rest (we later expect convergence to that very evec)
targetIndex = randi(N);
weights = rand(N,1);
weights(targetIndex) = weights(targetIndex) + weight_target;
targetV = V(:,targetIndex);
targetE = D(targetIndex);

v = V*weights;
v = v / norm(v);
fprintf("Matrix size: %d\n\n", N);

disp(['Target eigenvalue: ', num2str(targetE)]);


disp(['Residual = ', num2str(norm(a*v - rayleighquotient(a,v)*v))]);

disp(['Angle between v and target (deg): ', num2str(rad2deg(acos(v'* targetV)))]);
disp(['Angle between v and target (rad): ', num2str((acos(v'* targetV)))]);
disp(' ');

% sqrt(1 - (v'*targetV)^2)

% gamma = 1;
% atilde = a - gamma*1i*(speye(N) - v*v');
% [Vtilde, Dtilde] = eigs(atilde, N);
% Dtilde = diag(Dtilde); % extract eigenvalues

% Test classic RQI
tic; [x1, e1, its1, e_its1, v_its1, res1] = classic_rqi(a, v, inf); toc;
disp('Result of classic RQI');
disp(['Computed eigenvalue: ', num2str(e1), ' (', num2str(its1), ' Iterations)']);
if abs(targetE - e1) < 10e-9
    disp("Classic RQI succeded");
else
    disp("Classic RQI did not succeed");
end

disp(' ');

% Test complex RQI with gamma = res
tic; [x2, e2, its2, e_its2, v_its2, res2] = complex_rqi2(a, v, inf, inf); toc;
disp('Result of complex RQI (gamma = res)');
disp(['Computed eigenvalue: ', num2str(e2), ' (', num2str(its2), ' Iterations)']);
if abs(targetE - e2) < 10e-9
    disp("Complex RQI succeded");
else
    disp("Complex RQI did not succeed");
end

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

% The iteration number of V3 (squared residual) might sometimes be
% very high (> 1000), so we exclude it from the plots.
maxits = 40;
includeV3 = true;
if its3 > maxits
    includeV3 = false;
end

% Plot residuals
semilogy(1:its1+1, res1, '-o'); hold on
%loglog(1:its2+1, res2, '-x')
%if includeV3
%    loglog(1:its3+1, res3, '-s')    
%end
semilogy(1:its4+1, res4, '-^')
pbaspect([1 1 1])

xlabel("Iteration"); ylabel("Residual norm");


if includeV3
    max_x = max([its1,its2,its3,its4]);
else
    max_x = max([its1,its2,its4]);
end

max_x = max([its1, its4]);

min_y = min([res1, res2, res3, res4])/100;
max_y = max([res1, res2, res3 ,res4])*100;

axis([0, max_x+2, min_y, max_y]);

legend({'Classic RQI', 'Complex RQI'});

% if includeV3
%     legend({'Classic RQI', ...
%             'Complex RQI (\gamma^{(k)} = r)', ...
%             'Complex RQI (\gamma^{(k)} = r^2)',...
%             'Complex RQI (\gamma^{(k)} adaptive)'}, ...
%             'Location', 'southwest');
% else
%     legend({'Classic RQI', ...
%             'Complex RQI (\gamma^{(k)} = r)', ...
%             'Complex RQI (\gamma^{(k)} adaptive)'}, ...
%             'Location', 'southwest');
% end


%export_fig 'crqi_residuals.eps' -transparent