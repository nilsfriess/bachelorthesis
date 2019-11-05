clear;

disp('Close eigenvalues test');
disp(' ')

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

[testvecs,d] = eigs(Ahat, 20);
eigvals = diag(d);
eval10 = eigvals(10);
eval11 = eigvals(11);
eval12 = eigvals(12);

disp(['10th eigenvalue: e10 = ', num2str(eval10)])
disp(['11th eigenvalue: e11 = ', num2str(eval11)]) 
disp('Goal: Convergence to e10 and associated evec v10');
sigma = eval11 + 0.1;

disp(' ')
disp('TEST #1 (Shift left of e11 but initial vec close to v10')
disp(['Using classic RQI with shift s = ', num2str(sigma)])

evec10 = testvecs(:,10);
evec11 = testvecs(:,11);

n = 20;
t1 = linspace(0.05,0.95,n);
t2 = 1 - t1;

[r,~] = size(evec10);
testvecs = zeros(r,n);
res1 = zeros(n,4);

testvecs = normrnd(0,1,[r,n]);
testvecs = testvecs * nthroot(rand, n) / norm(testvecs);

%testvecs = testvecs';

for k = 1 : n
    %testvecs(:,k) = t1(k) * evec10 + (1 - t1(k)) * evec11;
    %testvecs(:,k) = testvecs(:,k) / norm(testvecs(:,k));
    % disp(['Angle between vtest and v10 in deg = ', num2str(rad2deg(acos(evec10'*testvecs(:,k))))])
    % disp(['Angle between vtest and v11 in deg = ', num2str(rad2deg(acos(evec11'*testvecs(:,k))))])
    angle = rad2deg(acos(evec10'*testvecs(:,k)));
    [vv, ee] = classic_rqi(Ahat, testvecs(:,k), sigma, 10e-8);
    [vc, ec] = complex_rqi(Ahat, testvecs(:,k), sigma, 2, 10e-8);
    % disp(['Angle to exact evec: ', num2str(angle), '  Converged to eigenvalue: ', num2str(ee)]);
    
    res1(k,1) = angle;
    res1(k,2) = ee;
    res1(k,3) = ec;
end

res1(:,4) = eval10;

sortrows(res1)
hold off
plot(res1(:,1), res1(:,2), 'b')
hold on
plot(res1(:,1), res1(:,3), 'r')
plot(res1(:,1), res1(:,4), 'y');



% [vv, ee] = classic_rqi(Ahat, vtest, sigma, 10e-8)