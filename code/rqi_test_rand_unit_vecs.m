% create 7x7 symmetric matrix with close eigenvalues
ldiag = rand([1,6]);
a = diag([4;12;15;41;40;60;90]) + diag(rand([1,7])) + diag(ldiag,1) + diag(ldiag,-1);

[V,D] = eigs(a, 7);
D = diag(D);

exact_eigval = D(3)
exact_eigvec = V(:,3)
exact_eigvec(3,1);

norm(a*exact_eigvec - exact_eigval*exact_eigvec)

% create n random unit vectors (will be the initial eigvec guesses)
n = 1000;
v = create_random_unit_vecs(7, n);

eigval_guess = exact_eigval + 0.5;
E = [];
V = [];
W = [];
V3 = [];

E_C = [];
V_C = [];
W_C = [];
V3_C = [];

for k = 1:n
    [x, sigma, its, eigval_its, eigvec_its] = classic_rqi(a, v(:,k), 40);
    V = [V, x];
    E = [E, sigma];
    W = [W, exact_eigvec'*v(:,k)];
    V3 = [V3, v(3,k)];
    
    [x, sigma, its, eigval_its, eigvec_its] = complex_rqi(a, v(:,k), 40, inf);
    V_C = [V_C, x];
    E_C = [E_C, sigma];
    W_C = [W_C, exact_eigvec'*v(:,k)];
    V3_C = [V3_C, v(3,k)];
end

DATA = [rad2deg(acos(W')), V3', E', E_C'];
DATA = sortrows(DATA, 1);

txt = sprintf("Eigenvalue to compute %f", exact_eigval);
disp(txt);

VarNames = { 'angle', 'v3', 'ClassicRQI', 'ComplexRQI' };
table(DATA(:,1), DATA(:,2), DATA(:,3), DATA(:,4), 'VariableNames', VarNames)
