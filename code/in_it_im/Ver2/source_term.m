function y = source_term(x1,x2)

y = (3*x1+x1^2)*exp(x1)*x2*(1-x2) + 2*x1*(1-x1)*exp(x1);