disp('Close eigenvalues test');
clear;
TEST = input("Choose Test: ");
disp(' ')

topright = 1.1;
[nodes,elements] = initial_mesh(topright);

for i = 1:5
    [nodes,elements] = refine(nodes,elements);
end

% First assemble the element stiffness matrices
elt_matrices = elt_stiffness(elements, nodes);

% Next assemble the global stiffness matrix
% (including all boundary nodes)
Ahat = global_stiffness(elt_matrices, elements, nodes);
% Ahat = load("../refine6.mat");
% Ahat = Ahat.a;
[rows,~] = size(Ahat);
disp(["Size of A" , rows]);

% Compute full eigendecomposition
[V,d] = eigs(Ahat, rows);
eigvals = diag(d);
eval10 = eigvals(10);
eval11 = eigvals(11);
eval12 = eigvals(12);
evec10 = V(:,10);
evec11 = V(:,11);

if TEST == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #1 (Shift on wrong side of wrong eigenvalue e11')
    disp('         and initial vector linear combination')
    disp('         of evec10 and evec11')
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp(['11th eigenvalue: e11 = ', num2str(eval11)])
    disp('Goal: Convergence to e10 and associated evec v10');
    
    sigma = eigvals(floor(1 + (rows - 1) * rand(1,1))) - 0.01;
    sigma = inf;
    lambda = 20;
    
    disp(['Using classic RQI with shift s = ', num2str(sigma), ' and'])
    disp(['Complex RQI with same real shift and complex shift ', num2str(lambda)])
    
    n = 100; % Number of tests
    t1 = linspace(0.01,0.99,n);
    t2 = 1 - t1;
    
    [r,~] = size(evec10);
    res1 = zeros(n,4);
    
    testvecs = zeros(r,n);
    weights = 0.2*rand(r,1);

    for k = 1 : n
        testvecs(:,k) = t1(k) * evec10 + (1 - t1(k)) * V(:,5);        
        testvecs = testvecs + (V'*weights);

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
    scatter(res1(:,1), res1(:,3), 'rx')
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
    sigma = eval11 - 0.1;
    lambda = 2;
    
    disp(['Using classic RQI with shift s = ', num2str(sigma), ' and'])
    disp(['Complex RQI with same real shift and complex shift ', num2str(lambda)])
    
    evec10 = V(:,10);
    n = 5;
    
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
    sigma = eval11 - 0.1;
    
    evec10 = V(:,10);
    n = 1000;
    
    lambda = linspace(0.1, 20, n);
    
    [r,~] = size(evec10);
    res1 = zeros(n,2);
    
    testvec = randn(r,1) + 6 * evec10;
    testvec = testvec / norm(testvec);
        
    [~, e] = classic_rqi(Ahat, testvec, sigma, 10e-8);
    
    for k = 1 : n
        [vc, ec] = complex_rqi(Ahat, testvec, sigma, lambda(k), 10e-12);
        
        res1(k,1) = lambda(k);
        res1(k,2) = ec;
    end
    res1(:,3) = eval10;
    
    disp(['Angle = ', num2str(rad2deg(acos(evec10'*testvec)))]);

    % sortrows(res1);
    figure
    plot([0, 20], [eval10, eval10], 'y');
    hold on
    plot([0, 20], [e, e], 'g');
    scatter(res1(:,1), res1(:,2), 'r')
    title('Plot of imag shift and convergence result');
    xlabel('shift');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}','RQI', 'CRQI'},'Location','best')
    
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
    sigma = eval11 - 0.01;
    lambda = 5;
    
    disp(['Using classic RQI with shift s = ', num2str(sigma), ' and'])
    disp(['Complex RQI with same real shift and complex shift ', num2str(lambda)])
    
    evec10 = V(:,10);
    
    weights = rand(n,rows);
    vecrand = V * weights(1,:)';   
    
        
    res1 = zeros(n,3);
    
    t1 = linspace(0.01,0.99,n);
    t2 = 2 - t1;
    
            
    for k = 1 : n
        
        testvec = t1(k) * vecrand + t2(k) * evec10;
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
    
elseif TEST == 5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #5 (Shift random between smallest and largest eigval')
    disp('         and initial vector random linear combination')
    disp('         of all eigenvectors, strongest in direction v10')
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp(['11th eigenvalue: e11 = ', num2str(eval11)])
    disp('Goal: Convergence to e10 and associated evec v10');
    
    lambda = 2;
        
    evec10 = V(:,10);
    
    n = 100; % Number of tests
    
    res1 = zeros(n,4);
    
    fprintf('sigma\t\t angle\t\t rqi\t\t crqi\t\t itsclassic\t itscomplex\n');
    
    for k = 1 : n
        weights = rand(rows,1);
        testvec = V * weights;
        testvec = testvec + 2.5 * evec10;
        testvec = testvec / norm(testvec);
                
        sigma = eigvals(1) * rand(1,1);

        anglev10 = rad2deg(acos(evec10'*testvec));
        anglev20 = rad2deg(acos(evec11'*testvec));
        [vv, ee, itsclassic] = classic_rqi(Ahat, testvec, inf, 10e-12);
        [vc, ec, itscomplex] = complex_rqi(Ahat, testvec, inf, lambda, 10e-12);
        
        res1(k,1) = sigma;
        res1(k,2) = ee;
        res1(k,3) = ec;
        
        fprintf("%f\t %f\t %f\t %f\t %f\t %f\n", sigma, anglev10, ee, ec, itsclassic, itscomplex);
    end
    res1(:,4) = eval10;
    
    % sortrows(res1);
    figure
    plot([min(res1(:,1)), max(res1(:,1))], [eval10, eval10], 'y');
    hold on
    scatter(res1(:,1), res1(:,2), 'bx')
    scatter(res1(:,1), res1(:,3), 'r')
    title('Plot of sigma and convergence result');
    xlabel('sigma');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}', 'RQI','CRQI'},'Location','best')
    
elseif TEST == 6
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #6 (Shift on wrong side of wrong eigenvalue e11')
    disp('         and initial vector random linear combination')
    disp('         of all eigenvectors, strongest in direction v10')
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp(['11th eigenvalue: e11 = ', num2str(eval11)])
    disp('Goal: Convergence to e10 and associated evec v10');
    
    sigma = 2.51; 
    lambda = 2;
        
    evec10 = V(:,10);
    
    n = 200; % Number of tests
    
    res1 = zeros(n,4);
    
    out = zeros(n,7);
    
    for k = 1 : n
        weights = rand(rows,1);
        testvec = V * weights;
        testvec = testvec + 1 * evec10;
        testvec = testvec / norm(testvec);
        
        anglev10 = rad2deg(acos(evec10'*testvec));
        anglev20 = rad2deg(acos(evec11'*testvec));
        [vv, ee] = classic_rqi(Ahat, testvec, sigma, 10e-8);
        [vc, ec] = complex_rqi(Ahat, testvec, sigma, lambda, 10e-8);
        
        res1(k,1) = anglev10;
        res1(k,2) = ee;
        res1(k,3) = ec;
        
        angles = sort(rad2deg(acos(V' * testvec)));
        
        
        resid = norm((Ahat - rayleighquotient(Ahat, testvec)*eye(rows))*testvec);
        gap = (eval10 - eval11)/2;
        out(k,:) = [angles(2), angles(2) - angles(1), anglev10, gap, resid, ee, ec];        
    end
    disp(num2str(sortrows(out, 2)))
    
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
    
elseif TEST == 7
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TEST #7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp(' ')
    disp('TEST #7 (Random real shift, fixed imag shift, fixed initvec)');
    
    disp(['10th eigenvalue: e10 = ', num2str(eval10)])
    disp(['11th eigenvalue: e11 = ', num2str(eval11)])
    disp('Goal: Convergence to e10 and associated evec v10');
    
    lambda = inf;
    disp(['Using complex shift lambda = ', num2str(lambda)]);
        
    angle = 100;
    max_search_iterations = 1000;
    while angle > 45
        evec10 = V(:,10);
        weights = rand(rows,1);
        testvec = V*weights;
        testvec = testvec + 4 * evec10;
        testvec = testvec / norm(testvec);
        angle = rad2deg(acos(evec10'*testvec));
        
        max_search_iterations = max_search_iterations - 1;
        if max_search_iterations == 0
           error("Did not find suitable starting vector");
           return;
        end
    end
    
    evalfirst = eigvals(1);
    evallast  = eigvals(end);
    
    n = 200; % Number of tests
    
    res1 = zeros(n,4);
    
    out = zeros(n,9);
    
    for k = 1 : n        
        anglev10 = acos(evec10'*testvec);
        
        sigma = evalfirst * rand(1,1);
        
        [vv, ee, itsr] = classic_rqi(Ahat, testvec, sigma, 10e-9);
        [vc, ec, itsc] = complex_rqi(Ahat, testvec, sigma, lambda, 10e-9);
        
        res1(k,1) = sigma;
        res1(k,2) = ee;
        res1(k,3) = ec;
        
        angles = sort(acos(V' * testvec));
        
        resid = norm((Ahat - sigma*eye(rows))*testvec);
        out(k,:) = [angles(2), angles(2) - angles(1), anglev10, sigma, resid, ee, ec, itsr, itsc];        
    end
    disp(num2str(sortrows(out, 4)))
    
    res1(:,4) = eval10;
    
    % sortrows(res1);
    figure
    plot([0, max(res1(:,1))], [eval10, eval10], 'y');
    hold on
    scatter(res1(:,1), res1(:,2), 'bx')
    scatter(res1(:,1), res1(:,3), 'r')
    title('Plot of initial shift and convergence result');
    xlabel('Initial shift');
    ylabel('Eigenvalue');
    legend({'\lambda_{10}', 'RQI','CRQI'},'Location','best')

    
end