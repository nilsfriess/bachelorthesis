function e = eigval_of_discrete_laplace(m, k)
%% Returns the kth eigenvalue of the matrix corresponding
%% to a discretized Laplacian. If k is a vector of integers, 
%% returns the eigenvalues given in the vector

    e = 4 * m^2 * (sin( (k.*pi)/(2*m) ) ).^2;
    e = e';
end