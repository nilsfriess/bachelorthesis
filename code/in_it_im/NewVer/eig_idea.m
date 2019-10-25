
function eig_idea(topright,nref)

% Arguments:
%             topright ... y coordinate of top-right corner (~ 1.0)
%             nref     ... number of refinements of initial mesh (< 6)
%             ilam     ... target eigenvalue (4, 5, or 6)
%             scale    ... scaling factor for initial guess of eval

format long
ii = j;

[nodes,elements] = initial_mesh(topright);
disp('Initial mesh is displayed in Figure Window')
visualise(nodes,elements)
disp('Type Return to Continue')
pause

disp('Refinement Level...')
for i = 1:nref
    [nodes,elements] = refine(nodes,elements);
    disp(i)
end
[nnodes,m] = size(nodes);

disp('Refined mesh is displayed in Figure Window (if not too large)')
if (nnodes < 1000)
    visualise(nodes,elements)
end
disp('Type Return to Continue')
pause

% First assemble the element stiffness matrices

elt_matrices = elt_stiffness(elements,nodes);

disp('Element stiffness matrices computed')

% Next assemble the global stiffness matrix 
% (including all boundary nodes)

Ahat = global_stiffness(elt_matrices,elements,nodes);

disp('Global stiffness matrix computed')

[nodes,index] = sortrows(nodes);

Ahat = Ahat(index,index);

N0 =  find(nodes(:,3)~=1);  % Finds the indices of the  
                            % interior (non-Dirichlet) nodes

A = Ahat(N0,N0);            % Select the appropriate blocks of Ahat
                            % to get the matrix A with Dirichlet BCs
                        
[V,D] = eig(full(A));

% Take interior eigenvalue lam_ilam as target eigenvalue and match it 
% with the second eigenvector  

N = length(N0)
    
% Set initial guess for target eigenvector

u = zeros(N,1);

Npos = find((nodes(N0,1)-0.5).^2+(nodes(N0,2)-0.5).^2<0.001);
   
u(Npos) = 1;
u = u/norm(u);

% Visualise initial guess

Uhat = zeros(nnodes,1);  
Uhat(N0) = u;

% Interpolate onto uniform mesh and plot initial guess using surf

h_unif = 1/(2^(nref+1));

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

surf(XI,YI,ZI);
disp('Type Return to Continue')
pause;

% Compute components of initial guess in direction of eigenvectors and
% display the largest ones

[Y,ind]=sort(abs(V'*u),1,'descend');

component = ind(1:5)
contrib = Y(1:5)

ilam = ind(1)
lamex = D(ilam,ilam)
vex = V(:,ilam);
    
% Display spectrum around target eigenvalue

spectrum = diag(D(max(1,ilam-3):ilam+3,max(1,ilam-3):ilam+3))

% Visualise target eigenvector

Uhat = zeros(nnodes,1);   % Define an extended vector 
Uhat(N0) = vex;

% Interpolate onto uniform mesh and plot exact solution using surf

h_unif = 1/(2^(nref+1));

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

surf(XI,YI,ZI);
disp('Type Return to Continue')
pause;

% Set initial guess for target eigenvalue (scaled by factor "scale")

%lam = scale*lamex

% Often initial guess for eigenvalue would be Rayleigh quotient

lam = u'*A*u

% Show how bad this guess is

plot(diag(D),zeros(N),'b.');
hold on
plot(lamex,0,'r*');
plot(lam,0,'g*');
hold off

% Compute initial residual and initial eval and evec errors

i=0
res = norm(A*u-lam*u)

lamerr = abs(lam-lamex)
uerr = 1 - abs(u'*vex)

disp('Type Return to Continue')
pause;

for i=1:100

    sig = lam - res^2*ii;
%    sig = lam;

    B = A - sig*eye(N);

%    Alternatively we can add on an imaginary projection into the 
%    orthogonal complement of the current evec guess. However, these
%    two approaches are equivalent in exact arithmetic. Nevertheless
%    it is useful in order to show why the method works (looking at 
%    the spectrum of the shifted system).
%
%    B = (A-lam*eye(N)) + res^2*ii*(eye(N)-u*u');
%
%    if (i==1)
%        [V2,D] = eig(B);
%
%        x=real(diag(D));
%        y=imag(diag(D));
%    
%        plot(x,y,'b.')
%        pause;
%    end
    
    % Carry out inverse iteration step, normalise and compute eigenvalue
    % approximation
    
    z = B\u;
    u = z/norm(z);
    lam = u'*A*u;
    
    % Compute residual and eval and evec errors
    
    i
    resold = res;
    res = norm(A*u-lam*u)
    
    lamerr = abs(lam-lamex)
    uerr = 1 - abs(u'*vex)

    % Compute components of current iterate in directions of eigenvectors
    % of A and display the largest ones
    
    [Y,ind]=sort(abs(V'*u),1,'descend');

    component = ind(1:5)
    contrib = Y(1:5)

    if max(res) < 1.0e-6
        break
    end
    
    disp('Type Return to Continue')
    pause;
         
end

u = real(u);
u = u/norm(u);

% Visualise final solution

Uhat = zeros(nnodes,1);  
Uhat(N0) = u;

% Interpolate onto uniform mesh and plot final solution using surf

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

surf(XI,YI,ZI);

end
