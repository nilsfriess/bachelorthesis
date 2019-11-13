%%% EIG SHIFT TEST

a = full(gallery('gcdmat', 4));
v1 = [1,1,1,1]';
v1 = v1 / norm(v1);

[V,D] = eigs(a);
diag(D)

v1'*V(:,1)

disp('     gamma     angle v(a,1), v(b.1)');
for i = 0:20
   gamma = i;
   b = a + gamma*1i*(eye(4) - v1*v1');
   
   [VV,DD] = eigs(b);
   diag(DD)
   out = [gamma, V(:,1)'*real(VV(:,1))];
   disp(out);
end