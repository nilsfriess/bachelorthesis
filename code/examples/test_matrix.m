function [a,m] = test_matrix(type, n)
    switch type
        case 1
            a = sprandsym(n,0.33);
        case 2
            % [1, 2, 1] matrix
            diagonals = [ones(1,n);2*ones(1,n);ones(1,n)]';
            a = spdiags(diagonals, [-1,0,1], n, n);
            
        case 3
            % Wilkinson's matrix W^{+} (note: n does not denote the size)
            p = n;
            n = 2*p + 1;
            main_diag = abs(p + 1 - (1:n));
            off_diag = ones(1,n);
            diagonals = [off_diag; main_diag; off_diag]';
            a = spdiags(diagonals, [-1,0,1], n, n);
            
        case 4
            % Martin-Wilkinson matrix MW
            a = zeros(n,n);
            for i = 1:n
                for j = 1:n
                    if i == j
                        a(i,j) = 6;
                    end
                    if (i == j-1) || (j == i-1)
                        a(i,j) = -4;
                    end
                    if (i == j-2) || (j == i-2)
                        a(i,j) = 1;
                    end
                end
            end
            a(1,1) = 5; a(end,end) = 5;
            a = sparse(a);
            
        case 5
            % Laplace matrix. Note that n does not denote the size
            m = n;
            n = m^2;
            I = speye(m);
            diagonals = [-ones(1,m);4*ones(1,m);-ones(1,m)]';
            T = spdiags(diagonals, [-1,0,1], m, m);
            a = kron(I,T);
            
            I = ones(n,1);
            offI = spdiags([-I -I], [-m,m], n,n);
            a = a + offI;
    end
    [m,~] = size(a);
end

