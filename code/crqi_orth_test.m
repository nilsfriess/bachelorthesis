function crqi_orth_test

n = 100;

a = gallery('fiedler', 1:n);

% initial vector guess
v = rand(n,1);
v = v / norm(v);

v_b = v;
v_c = v;

% initial eval guess
sigma = 10;
e_b = sigma;
e_c = sigma;

res_b = norm((a - sigma*eye(n))*v_b);
res_c = norm((a - sigma*eye(n))*v_c);

iteration = 0;


disp('  iter    e_b         e_c         r_b         r_c');

while (res_b > 10e-6) || (res_c > 10e-6)
    lambda_b = e_b - res_b^2*1i;
    lambda_c = e_c - res_c^2*1i;
    % assemble shifted systems
    b = a - lambda_b*(eye(n) - v_b*v_b');
    c = a - lambda_c*eye(n);
    
    v_b = b \ v_b;
    v_c = c \ v_c;
    
    v_b = v_b / norm(v_b);
    v_c = v_c / norm(v_c);
    
    e_b = rayleighquotient(a, v_b);
    e_c = rayleighquotient(a, v_c);
    
    res_b = norm((a - real(e_b)*eye(n))*real(v_b));
    res_c = norm((a - real(e_c)*eye(n))*real(v_c));
    
    out = [ iteration, e_b, e_c, res_b, res_c];
    disp(num2str(real(out)));
        
    iteration = iteration + 1;
end
    
end

