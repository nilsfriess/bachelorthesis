function E = laplace_eigvals(m)
    i = 1;
    E = zeros(m^2,1);
    for k=1:m
        for l=1:
            E(i) = (k^2 + ((2*l - 1)^2)/4)*pi^2;
            i = i + 1;
        end
    end
    E = sort(E);
end

