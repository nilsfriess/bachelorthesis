% program fem0.m  
%
% This solves the PDE problem (D) with homogeneous Dirichlet
% conditions using piecewise linear finite elements.

clear
format long

[nodes,elements] = initial_mesh(1.1);
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

N0 =  find(nodes(:,3)~=1);  % Finds the indices of the  
                            % interior (non-Dirichlet) nodes

A  = Ahat(N0,N0);       % Select the appropriate blocks of Ahat
                        % to get the matrix A from lectures

% Assemble the right hand side vector (the function in the RHS of
% (D) has to be specified in source_term)

f = rhs('source_term',N0,elements,nodes); 

U = A\f;  % Find coefficient vector U of finite element solution
          % at the degrees of freedom

% Finally for the purpose of drawing graphs of the solution
% and for computing errors, it is useful to assemble an extended  
% vector Uhat which contains the elements of U at the degrees of 
% freedom together with the Dirichlet boundary condition at the 
% boundary nodes. The ordering of the elements of Uhat coincides 
% with the ordering of the nodes.

Uhat = zeros(nnodes,1);   % Define an extended vector 

Uhat(N0) = U;  % put in the elements of U at the degrees of freedom

% Plot the solution using the MATLAB function surf
% It is first necessary to interpolate Uhat on a uniform mesh.
% We use the MATLAB functions meshgrid and griddata to do this.

h_unif = 1/(2^(nref+1));

[XI,YI] = meshgrid([min(nodes(:,1)):h_unif:max(nodes(:,1))],...
                   [min(nodes(:,2)):h_unif:max(nodes(:,2))]);
ZI      = griddata(nodes(:,1),nodes(:,2),Uhat,XI,YI,'linear');

hold on
surf(XI,YI,ZI);
hold off

disp('Numerical solution is displayed in Figure Window')
disp('Type Return to Continue')
pause
             
% Evaluate the exact solution at the nodes

x1 = nodes(:,1);
x2 = nodes(:,2);

Uex = x1.*(1-x1).*exp(x1).*x2.*(1-x2);

% Calculate error at nodes and plot it

err = abs(Uhat - Uex);

ZI = griddata(nodes(:,1),nodes(:,2),err,XI,YI,'linear');

surf(XI,YI,ZI);

disp('Error is displayed in Figure Window')
disp('Type Return to Continue')
pause

% Calculate a measure for the error (i.e. the energy norm of the difference 
% between the FE solution u_h and the interpolant of the exact solution u) 

anorm_err = sqrt((Uhat-Uex)'*Ahat*(Uhat-Uex))
