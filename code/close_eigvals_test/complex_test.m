TEST = 2;

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
[rows,~] = size(Ahat);

% Compute full eigendecomposition
[V,d] = eigs(Ahat, rows);
eigvals = diag(d);
eval10 = eigvals(10);
eval11 = eigvals(11);
eval12 = eigvals(12);

if TEST == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #1 (Shift right of e11 (wrong eigval)')
    disp('         and initial vector linear combination')
    disp('         of evec10 and evec11')
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp(['11th eigenvalue: e11 = ', num2str(eval11)])
    disp('Goal: Convergence to e10 and associated evec v10)');
    sigma = eval11-0.001;
    lambda = 5;
    
    disp(['Using classic RQI with shift s = ', num2str(sigma), ' and'])
    disp(['Complex RQI with same real shift and complex shift ', num2str(lambda)])
    
    evec10 = V(:,10);
    evec11 = V(:,11);
    
    n = 1000;
    t1 = linspace(0.01,0.99,n);
    t2 = 1 - t1;
    
    [r,~] = size(evec10);
    res1 = zeros(n,4);
    
    testvecs = zeros(r,n);
    
    for k = 1 : n
        testvecs(:,k) = t1(k) * evec10 + (1 - t1(k)) * evec11;
        testvecs(:,k) = testvecs(:,k) / norm(testvecs(:,k));
        testvecs(:,k) = testvecs(:,k) / norm(testvecs(:,k));
        angle = rad2deg(acos(evec10'*testvecs(:,k)));
        [vv, ee] = classic_rqi(Ahat, testvecs(:,k), sigma, 10e-8);
        [vc, ec] = complex_rqi(Ahat, testvecs(:,k), sigma, lambda, 10e-8);
        
        res1(k,1) = angle;
        res1(k,2) = ee;
        res1(k,3) = ec;
    end
    res1(:,4) = eval10;
    
    % sortrows(res1);
    figure
    plot([0, 90], [eval10, eval10], 'y');
    hold on
    plot(res1(:,1), res1(:,2), 'b')
    plot(res1(:,1), res1(:,3), 'r')
    title('Plot of angle(v_{10}, initv) and convergence result');
    xlabel('angle');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}', 'RQI','CRQI'},'Location','best')
    
elseif TEST == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #2 (Shift right of e11 (wrong eigval)')
    disp('         and initial vector linear combination')
    disp('         of evec10 and random unit vec')
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp('Goal: Convergence to e10 and associated evec v10)');
    sigma = eval11 - 0.001;
    lambda = 2;
    
    disp(['Using classic RQI with shift s = ', num2str(sigma), ' and'])
    disp(['Complex RQI with same real shift and complex shift ', num2str(lambda)])
    
    evec10 = V(:,10);
    n = 1000;
    
    [r,~] = size(evec10);
    res1 = zeros(n,4);
    
    testvecs = zeros(r,n);
    
    for k = 1 : n
        testvecs(:,k) = 0.2*randn(r,1); %+ 0.9*evec10;
        testvecs(:,k) = testvecs(:,k) / norm(testvecs(:,k));
        testvecs(:,k) = testvecs(:,k) / norm(testvecs(:,k));
        angle = rad2deg(acos(evec10'*testvecs(:,k)));
        [vv, ee] = classic_rqi(Ahat, testvecs(:,k), sigma, 10e-8);
        [vc, ec] = complex_rqi(Ahat, testvecs(:,k), sigma, lambda, 10e-8);
        
        res1(k,1) = angle;
        res1(k,2) = ee;
        res1(k,3) = ec;
    end
    res1(:,4) = eval10;
    
    % sortrows(res1);
    figure
    plot([0, 90], [eval10, eval10], 'y');
    hold on
    scatter(res1(:,1), res1(:,2), 'b')
    scatter(res1(:,1), res1(:,3), 'r')
    title('Plot of angle(v_{10}, initv) and convergence result');
    xlabel('angle');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}', 'RQI','CRQI'},'Location','best')
    
elseif TEST == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #3 (Shift right of e11 (wrong eigval)')
    disp('         initial vector fixed')
    disp('         and varying imaginary shift in CRQI')
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp('Goal: Convergence to e10 and associated evec v10)');
    sigma = eval11 - 0.01;
    
    evec10 = V(:,10);
    n = 1000;
    
    lambda = linspace(0.1, 20, n);
    
    [r,~] = size(evec10);
    res1 = zeros(n,2);
    
    testvec = randn(r,1) + evec10;
    testvec = testvec / norm(testvec);
    
    disp(['Angle = ', num2str(rad2deg(acos(evec10'*testvec)))]);
    
    [~, e] = classic_rqi(Ahat, testvecs(:,k), sigma, 10e-8);
    
    for k = 1 : n
        [vc, ec] = complex_rqi(Ahat, testvecs(:,k), sigma, lambda(k), 10e-8);
        
        res1(k,1) = lambda(k);
        res1(k,2) = ec;
    end
    res1(:,3) = eval10;
    
    % sortrows(res1);
    figure
    plot([0, 20], [eval10, eval10], 'y');
    hold on
    scatter(res1(:,1), res1(:,2), 'r')
    plot([0,20], [e,e], 'b');
    title('Plot of imag shift and convergence result');
    xlabel('shift');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}','CRQI', 'RQI'},'Location','best')
    
elseif TEST == 4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #4 (Shift right of e11 (wrong eigval)')
    disp('         and initial vector linear combination')
    disp('         of all eigenvectors, but strongest')
    disp('         in the direction of v10')
    
    n = 1000;
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp('Goal: Convergence to e10 and associated evec v10)');
    sigma = eval11-0.005;
    lambda = 10;
    
    disp(['Using classic RQI with shift s = ', num2str(sigma), ' and'])
    disp(['Complex RQI with same real shift and complex shift ', num2str(lambda)])
    
    evec10 = V(:,10);
    
    weights = rand(1,rows);
    weights10 = linspace(0.5,10,n);
        
    res1 = zeros(n,3);
    
            
    for k = 1 : n
        testvec = zeros(rows, 1);
        weights(10) = weights10(k);
        for j = 1 : rows
            testvec = testvec + V(:,j)*weights(j);
        end
        testvec = testvec / norm(testvec);
        
        angle = rad2deg(acos(evec10'*testvec));
        [vv, ee] = classic_rqi(Ahat, testvec, sigma, 10e-8);
        [vc, ec] = complex_rqi(Ahat, testvec, sigma, lambda, 10e-8);
        
        res1(k,1) = angle;
        res1(k,2) = ee;
        res1(k,3) = ec;
    end
    res1(:,4) = eval10;
    
    % sortrows(res1);
    figure
    plot([0, 90], [eval10, eval10], 'y');
    hold on
    scatter(res1(:,1), res1(:,2), 'b')
    scatter(res1(:,1), res1(:,3), 'r')
    title('Plot of angle(v_{10}, initv) and convergence result');
    xlabel('angle');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}', 'RQI','CRQI'},'Location','best')
end