function DrazinInverse1a = DrazinInverse(a)
%-----------------------------------------
%Compute the Drazin Inverse of a matrix 'a' using the limited algorithm.
%Need computing the index of 'a'.
global q1 q2 s1 s2
[m,n] = size(a);
if m~= n
    disp('Matrix is must be square!')
end
%-----------------------------------------
% Computer the index of A and note r = rank(A^k).
[k,r,a1,a] = index(a);
F = eye(n);
g = -trace(a);
%g = collect(g);
for i = 1:r-1
    G = g*eye(n);
    F = a*F+G;
    g = -trace(a*F)/(i+1)
    %g = collect(g);
end
DrazinInverse1a = a1*F;
DrazinInverse1a = -1/g*DrazinInverse1a;
%DrazinInverse1a = simplify(DrazinInverse1a);
end

function [k,r,a1,a] = index(a)
%To compute the index of 'a'.
k = 0;
n = length(a);
r = n;
a0 = a;
r1 = rank(a);
a1 = eye(n);
while r ~= r1
    r = r1;
    a1 = a;
    a = a*a0;
    r1 = rank(a);
    k = k+1;
end
r
%r = sym2poly(r);
end