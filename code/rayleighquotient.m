function r = rayleighquotient(a,x)
%RAYLEIGHQUOTIENT Computes the Rayleigh Quotient of x wrt to the matrix a
    x = x / norm(x);
    r = x'*a*x;
end

