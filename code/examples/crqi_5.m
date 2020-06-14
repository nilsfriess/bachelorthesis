clf;

test = 2;
size = 1000;
weight_target_max = 100;
 
y_d_min = 0.5;
y_d_max = 0.5;

tik_up = 0.22;


[a,N] = test_matrix(test,size);
fprintf("Matrix size: %d\n\n", N);

[V,D] = eigs(a,N);
D = diag(D);

tests = 250;
addWeights = linspace(0,weight_target_max,tests);

% Create random linear combination of eigenvectors of A
% where one of the interior eigenvectors is weighted heavier
% than the rest (we later expect convergence to that very evec)
%targetIndex = randi([1,N]);
targetIndex = 100;
targetV = V(:,targetIndex);
targetE = D(targetIndex);

disp(['Target eigenvalue: ', num2str(targetE)]);

resRQI = zeros(tests,3);
resCRQI= zeros(tests,3);

weights = rand(N,1);
weights(targetIndex) = weights(targetIndex) + addWeights(end);

v = V*weights;
v = v / norm(v);
fprintf("Minimum angle: %.5f\n\n", rad2deg(acos(abs(v'*targetV))));


for i = 1:tests
    fprintf('# %4d  ', i);
    
    weights = rand(N,1);
    weights(targetIndex) = weights(targetIndex) + addWeights(i);
    
    v = V*weights;
    v = v / norm(v);
    
    
    angle = acos(abs(v'* targetV));
    
    fprintf("Angle: %.5f, ", rad2deg(angle));
    
    [x1, e1, its1, e_its1, v_its1, res1] = classic_rqi(a, v, inf);
    [x2, e2, its2, e_its2, v_its2, res2] = complex_rqi2combined(a, v, inf, inf);
    
    if abs(e1 - targetE) < 10e-8
        resRQI(i,3) = true;
    else
        resRQI(i,3) = false;
    end
    
    if abs(e2 - targetE) < 10e-8
        resCRQI(i,3) = true;
    else
        resCRQI(i,3) = false;
    end
    
    resRQI(i,1) = angle;
    resCRQI(i,1) = angle;

    resRQI(i,2) = its1;
    resCRQI(i,2) = its2;
    
    fprintf("Iteration count: RQI = %2d, CRQI = %2d\n", its1, its2);
end

% Split into success and failure
rqiSuccess = resRQI(resRQI(:,3) == true, :);
rqiFail = resRQI(resRQI(:,3) == false, :);

crqiSuccess = resCRQI(resCRQI(:,3) == true, :);
crqiFail = resCRQI(resCRQI(:,3) == false, :);

% To be able to distinguish the plots, move the result for 
% classic RQI slightly up
rqiSuccess(:,2) = rqiSuccess(:,2) + tik_up;
rqiFail(:,2) = rqiFail(:,2) + tik_up;

plot(rqiSuccess(:,1), rqiSuccess(:,2), 'ok'); hold on;
plot(rqiFail(:,1), rqiFail(:,2), 'xk');

plot(crqiSuccess(:,1), crqiSuccess(:,2), 'or');
plot(crqiFail(:,1), crqiFail(:,2), 'xr');

if isempty(rqiSuccess(:,2))
    legend(["Classic RQI, wrong result",...
        "Complex RQI, correct result",...
        "Complex RQI, wrong result"], 'Location', 'northwest');
else
    legend(["Classic RQI, correct result",...
        "Classic RQI, wrong result",...
        "Complex RQI, correct result",...
        "Complex RQI, wrong result"], 'Location', 'northwest');
end

max_y = max([max(crqiSuccess(:,2)), max(crqiFail(:,2))]) + 1;
    
axis([0,pi/2+0.1,1,max_y]);

xticklabels({'\pi/12', '\pi/6', '\pi/4', '\pi/3', '5\pi/12', '\pi/2'})
xticks([pi/12, pi/6, pi/4, pi/3, 5*pi/12, pi/2])

ylabel('Iterations')
pbaspect([ 1 1 1 ])
xlabel('Angle between initial vector and target')

fprintf("Matrix size: %d\n\n", N);
disp(['Target eigenvalue: ', num2str(targetE)]);

box on;





