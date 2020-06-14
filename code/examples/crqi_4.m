% Wilkinson matrix
p1 = 10;   n1 = 2*p1+1;  W1 = test_matrix(3, p1);
p2 = 52;   n2 = 2*p2+1;  W2 = test_matrix(3, p2);
p3 = 262;  n3 = 2*p3+1;  W3 = test_matrix(3, p3);
p4 = 514;  n4 = 2*p4+1;  W4 = test_matrix(3, p4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TEST #1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("Test #1, n = %d\n\n", n1);

[V1,E1] = eigs(W1, n1);
E1 = diag(E1);

e_1_5 = E1(5);
e_1_6 = E1(6);
v_1_5 = V1(:,5);
v_1_6 = V1(:,6);

fprintf("Target eigenvalue: %.5e\n", e_1_5);
fprintf("Difference between e5 and e6: %.5e\n", abs(e_1_5 - e_1_6));

weights = rand(n1,1);
weights(5) = 2;
u = V1*weights;
u = u / norm(u);

disp(['Angle: ', num2str(rad2deg(acos(u'*v_1_5)))]);
disp(['Init RQ: ', num2str(u'*W1*u)]);

[x1, e1, its1, e_its1, x_its1, res1] = complex_rqi2combined(W1, u, inf, inf);

[x2, e2, its2, e_its2, x_its2, res2] = classic_rqi(W1, u, inf);

