function eig_idea_1d(m)

% Set m to be the nearest multiple of 4 minus 1

ii = j;

m = ceil(m/8)*8-1;

A = (2*eye(m) - diag(ones(m-1,1),1) - diag(ones(m-1,1),-1));

[V,D] = eig(A);

diag(D((m+1)/4-3:(m+1)/4+3,(m+1)/4-3:(m+1)/4+3))

lam_ex = D((m+1)/4,(m+1)/4);
v_ex = V(:,2);

D((m+1)/4,(m+1)/4) = D(2,2);
D(2,2) = lam_ex;

A = V*D*V';

%[V,D] = eig(A);

% Start with a shift that is 50% wrong and an approximation of the
% eigenvector that has roughly the right features 

lam = 0.7*lam_ex

u = ones(m,1);
u((m+1)/2-m/4+1:(m+1)/2-m/8) = 0;
u(1:(m+1)/2-m/4) = -1;
u = u/norm(u);

plot(u);
hold on
plot(v_ex,'m-');
hold off
pause

lamerr = abs(lam-lam_ex)
res = norm(A*u-lam*u)

for i=1:20

    B = (A-lam*eye(m)) + res*ii*(eye(m)-u*u');
    
    [V,D] = eig(B);

    x=real(diag(D));
    y=imag(diag(D));
    
    plot(x,y,'b.')

    z = B\u;
    u = real(z/norm(z));
    lam = u'*A*u;
    
    i
    resold = res;
    res = norm(A*u-lam*u)
    res/resold^2
    
    lamerr = abs(lam-lam_ex);
    uerr = 1-abs(u'*v_ex);
    
    if max(res) < 1.0e-9
        break
    end
    
    pause;
         
end 

end

