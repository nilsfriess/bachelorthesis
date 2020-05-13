clf;

test = 5;
size = 100;
weight_target_max = 180;

y_d_min = 0.5;
y_d_max = 0.5;


[a,N] = test_matrix(test,size);
fprintf("Matrix size: %d\n\n", N);

% a = load("fe_matrix.mat");
% a = a.Ahat;
% [N,~] = size(a)

[V,D] = eigs(a,N);
D = diag(D);

tests = 250;
addWeights = linspace(0,weight_target_max,tests);

classicRQIResults = zeros(tests,5);
complexRQIResults = zeros(tests,5);

% Create random linear combination of eigenvectors of A
% where one of the interior eigenvectors is weighted heavier
% than the rest (we later expect convergence to that very evec)
targetIndex = randi(N);
targetV = V(:,targetIndex);
targetE = D(targetIndex);

disp(['Target eigenvalue: ', num2str(targetE)]);

for i = 1:tests
    disp(['Iteration: ', num2str(i)]);
    
    weights = rand(N,1);
    weights(targetIndex) = weights(targetIndex) + addWeights(i);
    
    v = V*weights;
    v = v / norm(v);
    
    Vwithouttarget = V;
    Vwithouttarget(:,targetIndex) = zeros(N,1);
    sumOfEvecs = sum(Vwithouttarget, 2);
    sumOfEvecs = sumOfEvecs / norm(sumOfEvecs);
    
    disp(['Angle between v and target (deg): ', num2str(rad2deg(acos(v'* targetV)))]);
    disp(['Angle between v and target (rad): ', num2str((acos(v'* targetV)))]);
    disp(['Angle between v and remain (deg): ', num2str(rad2deg(acos(v'* sumOfEvecs)))]);
    
    [x1, e1, its1, e_its1, v_its1, res1] = classic_rqi(a, v, inf);
    [x2, e2, its2, e_its2, v_its2, res2] = complex_rqi2combined(a, v, inf, inf);
    
    classicRQIResults(i,1) = acos(v'*targetV);    
    complexRQIResults(i,1) = acos(v'*targetV);
    
    classicRQIResults(i,2) = acos(v'*sumOfEvecs);    
    complexRQIResults(i,2) = acos(v'*sumOfEvecs);
    
    classicRQIResults(i,4) = its1;    
    complexRQIResults(i,4) = its2;
    
    classicRQIResults(i,5) = e1;    
    complexRQIResults(i,5) = e2;

    if abs(e1 - targetE) < 10e-4
        classicRQIResults(i,3) = 1;
    end
    
    if abs(e2 - targetE) < 10e-4
        complexRQIResults(i,3) = 1;
    end
end

% figure(1);
% subplot(2,2,1);
% plot(complexRQIResults(:,1), complexRQIResults(:,3), 'xk'); hold on;
% plot(classicRQIResults(:,1), classicRQIResults(:,3)+0.2, 'or');
% axis([min(complexRQIResults(:,1))-0.5,...
%       max(complexRQIResults(:,1))+0.5,...
%       -3,3]);
% 
% subplot(2,2,2);
% plot(complexRQIResults(:,2), complexRQIResults(:,3), 'xk'); hold on;
% plot(classicRQIResults(:,2), classicRQIResults(:,3)+0.2, 'or');
% axis([min(complexRQIResults(:,2))-0.5,...
%       max(complexRQIResults(:,2))+0.5,...
%       -3,3]);
% 
% subplot(2,2,3);
% plot(complexRQIResults(:,1), complexRQIResults(:,4), 'xk'); hold on;
% plot(classicRQIResults(:,1), classicRQIResults(:,4), 'or');
% 
% subplot(2,2,4);
% plot(complexRQIResults(:,2), complexRQIResults(:,4), 'xk'); hold on;
% plot(classicRQIResults(:,2), classicRQIResults(:,4), 'or');


plot(complexRQIResults(:,1), complexRQIResults(:,5), 'xk');
hold on;
plot(classicRQIResults(:,1), classicRQIResults(:,5), 'or');
plot([min(complexRQIResults(:,1)) - 1, pi/2 + 1], [targetE, targetE], 'k');
ylabel('Computed eigenvalue')
pbaspect([ 1 1 1 ])
xlabel('Angle between initial vector and target')



legend({'Complex RQI', 'Classic RQI', 'Target eigenvalue'})


axis([min(complexRQIResults(:,1)) - 0.1, pi/2 + 0.1,...
      min(complexRQIResults(:,5))-y_d_min,max(complexRQIResults(:,5))+y_d_max]);

xticklabels({'\pi/12', '\pi/6', '\pi/4', '\pi/3', '5\pi/12', '\pi/2'})
xticks([pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2])

fprintf("Matrix size: %d\n\n", N);
