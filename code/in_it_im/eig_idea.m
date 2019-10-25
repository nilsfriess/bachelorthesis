% program eig_idea 

clear
format long
ii = j;

topright = input('Choose top right corner of domain (something close to 1.0):');

[nodes,elements] = initial_mesh(topright);
disp('Initial mesh is displayed in Figure Window')
visualise(nodes,elements)
disp('Type Return to Continue')
pause

nref = input('Type how many times you want to refine the initial mesh:  ');

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

Aorig  = Ahat(N0,N0);       % Select the appropriate blocks of Ahat
                            % to get the matrix A from lectures
                        
[V,D] = eig(full(Aorig));

% Take interior eigenvalue lam_ilam as target eigenvalue and match it 
% with the second eigenvector  

N = length(N0)

ilam = 4; %input('Please choose index of target eigenvector')
%ilam = floor(N/6)
lamex = D(ilam,ilam)
vex = V(:,4);
%vex = V(:,1);
     
% Give spectrum around target eigenvalue (see double eigenvalues)

diag(D(ilam-3:ilam+3,ilam-3:ilam+3))
%diag(D(ilam:ilam+3,ilam:ilam+3))

% Visualise 2nd (target) eigenvector

Uhat = zeros(nnodes,1);   % Define an extended vector 
Uhat(N0) = vex;

% Interpolate onto uniform mesh and plot solution using surf

h_unif = 1/(2^(nref+1));

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

surf(XI,YI,ZI);
pause;

%D(ilam,ilam) = D(4,4);
%D(4,4) = lamex;

%A = V*D*V';
A = Aorig;

%spy(A>0.0001);

% Set initial guess for target eigenvalue

scale = input('Scale exact eigenvalue by what factor?');
lam = scale*lamex

% Show how bad this guess is

plot(diag(D),zeros(N),'b.');
hold on
plot(lamex,0,'r*');
plot(lam,0,'g*');
hold off
pause;

% Set initial guess for target eigenvector

u = zeros(N,1);

Npos = find((nodes(N0,1)<0.4 & nodes(N0,2)<0.5) | ...
            (nodes(N0,1)>0.6 & nodes(N0,2)>0.45) );
Nneg = find((nodes(N0,1)>0.45 & nodes(N0,2)<0.4)| ...
            (nodes(N0,1)<0.4 & nodes(N0,2)>0.65));

u(Npos) = 1;
u(Nneg) = -1;
%u = ones(N,1);
u = u/norm(u);

[Y,ind]=sort(abs(V'*u),1,'descend');

ind(1:10)
Y(1:10)

% lam = u'*A*u

% Visualise initial guess

Uhat = zeros(nnodes,1);  
Uhat(N0) = u;

% Interpolate onto uniform mesh and plot solution using surf

h_unif = 1/(2^(nref+1));

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

surf(XI,YI,ZI);
pause;

i=0
res = norm(A*u-lam*u)

lamerr = abs(lam-lamex)
uerr = 1 - abs(u'*vex)

pause;
sig = lam;
for i=1:200

    B = (A-sig*eye(N)) + res^2*ii*(eye(N)-u*u');

    if (i==1)
        [V,D] = eig(B);

        x=real(diag(D));
        y=imag(diag(D));
    
        plot(x,y,'b.')
    end
    
    z = B\u;
    %[z,flag,relres,iter] = gmres(B,u,[], min(res,0.001),N); %,[],[],u);
    %its = iter(2)
    
    %z = real(z);
    u = z/norm(z);
    lam = u'*A*u;
    
    i
    resold = res;
    res = norm(A*u-lam*u)
    convl(i)=res/resold;
    convq(i)=res/resold^2;
    convc(i)=res/resold^3;
    
    lamerr = abs(lam-lamex)
    uerr = 1 - abs(u'*vex)
    %uerr = acos(abs(u'*vex))
    
    if max(res) < 1.0e-6
        break
    end
    
    pause;
         
end
u = real(u);
u = u/norm(u);

plot(convl(1:i-2),'k*');
pause
plot(convq(1:i-2),'k*');
pause
plot(convc(1:i-2),'k*');
pause

% Visualise final solution

Uhat = zeros(nnodes,1);  
Uhat(N0) = u;

% Interpolate onto uniform mesh and plot solution using surf

h_unif = 1/(2^(nref+1));

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

surf(XI,YI,ZI);




