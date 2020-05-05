function a = test_matrix(type, n)
    switch type
        case 1
            % [1, 2, 1] matrix
            diagonals = [ones(1,n);2*ones(1,n);ones(1,n)]';
            a = spdiags(diagonals, [-1,0,1], n, n);
            
        case 2
            % Wilkinson's matrix W^{-} (note: n does not denote the size)
            p = n;
            n = 2*p + 1;
            main_diag = p + 1 - (1:n);
            off_diag = ones(1,n);
            diagonals = [off_diag; main_diag; off_diag]';
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
            I = eye(m);
            T = 4*diag(ones(m,1))-(diag(ones(m-1,1),1) + diag(ones(m-1,1),-1));
            a = kron(I,T);
                        
            offdiag = -1*ones(m^2 - m,1);
            for i = 1:m
                offdiag(end-i+1) = -sqrt(2);
            end
            a = a + (diag(offdiag, m) + diag(offdiag,-m));
            a = sparse(a);
    end
end

