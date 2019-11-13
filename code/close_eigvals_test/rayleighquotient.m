function r = rayleighquotient(a,x)
%RAYLEIGHQUOTIENT Computes the Rayleigh quotient of x wrt. the matrix a
    x = x / norm(x);
    r = x' * a * x;
end

