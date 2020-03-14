m = 100;
density = 0.4;
rc = 0.2;
a = sprandsym(100, 0.6, rc, 1);

clear('targetVec');
clear('otherEigVecs');

[V,D] = eigs(a, m);

col = 20;
% sum of all but one eigenvectors
otherEigVecs = sumOfColsWithout(a, col);
otherEigVecs = otherEigVecs / norm(otherEigVecs);

targetVec = V(:,col);

vals = 0:0.01:1;
[~, valSize] = size(vals);
res = zeros(valSize, 5);
res(:,4) = D(col,col);
ctr = 1;
for i = vals
    % test vector
    initVec = i * targetVec + (1-i) * otherEigVecs;
    initVec = initVec / norm(initVec);
    
    disp(['Angle between target and initvec: ', num2str(targetVec'*initVec)]);
    disp(['Angle between rest of spectrum and initvec: ', num2str(initVec'*otherEigVecs)]); 
    
    [x,e] = complex_rqi_p(a,initVec,inf,100);
    [x,e1] = classic_rqi(a,initVec,inf);
    disp(['Target eigenvalue: ', num2str(D(col,col))]);
    disp(['Computed eigenvalue: ', num2str(e)]);
    disp(' ');
    
    res(ctr,1) = abs(targetVec'*initVec);
    res(ctr,2) = abs(initVec'*otherEigVecs);
    res(ctr,3) = e;
    res(ctr,5) = e1;
    ctr = ctr + 1;
end

subplot(1,2,1);
plot([min(res(:,2)),max(res(:,2))], [res(1,4), res(1,4)], 'r');
xlabel("Angle between target and init");
hold on;
plot(res(:,1), res(:,3), 'o');
plot(res(:,1), res(:,5), 'x');

subplot(1,2,2);
plot([min(res(:,2)),max(res(:,2))], [res(1,4), res(1,4)], 'r');
xlabel("Angle between sum of other evecs and init");
hold on;
plot(res(:,2), res(:,3), 'o');
plot(res(:,2), res(:,5), 'x');

