function v = eigvec_of_discrete_laplace(m, k)
%% Returns the normalised eigenvector corresponding to the kth
%% eigenvalue of the matrix corresponding to a discretized
%% Laplacian 

    v = ones(m-1, length(k));
    for n = 1:length(k)
        for j = 1:m-1
            v(j,n) = sin( (n*pi*j)/m );
        end
        v(:,n) = v(:,n) / norm(v(:,n));
    end
end