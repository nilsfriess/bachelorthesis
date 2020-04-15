clf

N = 200;
a = sprandsym(N,0.2);

% a = load("fe_matrix.mat");
% a = a.Ahat;
% [N,~] = size(a)

[V,D] = eigs(a,N);
D = diag(D); % D is diagonal matrix, extract eigenvalues and store back

% Create random linear combination of eigenvectors of A
% where one of the interior eigenvectors is weighted heavier
% than the rest (we later expect convergence to that very evec)
targetIndex = randi(N);
weights = rand(N,1);
weights(targetIndex) = 3;
targetV = V(:,targetIndex);
targetE = D(targetIndex);

Vwithouttarget = V;
Vwithouttarget(:,targetIndex) = zeros(N,1);
sumOfEvecs = sum(Vwithouttarget, 2);
sumOfEvecs = sumOfEvecs / norm(sumOfEvecs);

disp(['Target eigenvalue: ', num2str(targetE)]);


v = V*weights;
v = v / norm(v);


% Random unit vector
v = randn(N,1);
v = v / norm(v);

disp(['Residual = ', num2str(norm(a*v - rayleighquotient(a,v)*v))]);

disp(['Angle between v and target (deg): ', num2str(rad2deg(acos(v'* targetV)))]);
disp(['Angle between v and target (rad): ', num2str((acos(v'* targetV)))]);
disp(['Angle between v and remain (deg): ', num2str(rad2deg(acos(v'* sumOfEvecs)))]);
disp(' ');


maxits = 100;
shifts = linspace(min(D)-10, max(D)+10, maxits);
%shifts = linspace(-20,20,maxits);

result_crqi = zeros(maxits,1);
result_rqi = zeros(maxits,1);
for k = 1 : maxits
    currshift = shifts(k);
    
    [x1, e1, its1, e_its1, v_its1, res1] = classic_rqi(a, v, currshift);
    [x2, e2, its2, e_its2, v_its2, res2] = complex_rqi2(a, v, currshift, inf);
    
    result_crqi(k) = e2;
    result_rqi(k) = e1;
end

% Make dashed area between lower and upper bound of spectrum
% [X,Y] = hatch_coordinates([min(D),max(D)], ...
%                           [0,(min(result_rqi)-2)], ...
%                           0.5, 0.5);
% plot(X,Y,'Color',[.7 .7 .7],'linewidth',1,'LineStyle',':');
% hold on;
% [X,Y] = hatch_coordinates([min(D),max(D)], ...
%                           [0, (max(result_rqi)+2)], ...
%                           0.5, 0.5);
% plot(X,Y,'Color',[.7 .7 .7],'linewidth',1,'LineStyle',':');

P = patch('XData', [min(D), min(D), max(D), max(D)],...
          'YData', [min(result_rqi) - 2, ...
                    max(result_rqi) + 2, ...
                    max(result_rqi) + 2, ...
                    min(result_rqi) - 2],...
          'FaceColor', 'White',...
          'FaceAlpha', 0,...
          'EdgeColor', [1 1 1]);
hh1 = hatchfill(P, 'single', 45, 6, [0.2 0.2 0.2]);
hh1.Color = [.3 .3 .3 0.7];
hh1.LineStyle = ':';
hold on;

axis([min(shifts) inf (min(result_rqi)-1)  (max(result_rqi)+1)]);
pbaspect([1 1 1])

p1 = plot(shifts, result_rqi, 'or'); 
p2 = plot(shifts, result_crqi, 'xb');
%p3 = plot([shifts(1), shifts(end)], [targetE, targetE], 'k');
%plot([min(D), max(D)], [targetE, targetE], '^k');
legend([p1,p2], {'Classic RQI', 'Complex RQI'}, 'Location', 'northwest')
xlabel('Initial real shift') 
ylabel('Computed eigenvalue')
box on

%export_fig 'crqi_varyingshift__rand_initvec.eps' -transparent 