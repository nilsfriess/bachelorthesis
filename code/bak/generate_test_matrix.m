function a = generate_test_matrix(m)
    a = diag(2 * ones(1,m-1));
    a = a + diag(-1 * ones(1,m-2), 1);
    a = a + diag(-1 * ones(1,m-2),-1);
   
    a = m^2 * a;   
end