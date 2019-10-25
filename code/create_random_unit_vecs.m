function [v] = create_random_unit_vecs(dim, n)
    v = zeros(dim, n);
    for i = 1:n
        x = zeros(dim,1);
        while norm(x) < .0001
            for m = 1:dim
                x(m) = randn();
            end
        end
    
        x = x / norm(x); % normalise
        v(:,i) = x;
    end
end

