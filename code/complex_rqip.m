function [x, sigma, iterations, eigval_iterates, eigvec_iterates, residuals] = complex_rqip(a,x,sigma,tolerance)
% complex_rqi   Computes an eigenpair of a using the complex Rayleigh 
% quotient iteration
%
%   See CLASSIC_RQI for description of return values. The parameters a, x, 
%   sigma and tol are also described there.

  if nargin < 3
    error("Too few arguments");
  elseif nargin < 4
    tolerance = 10e-9;
  end
  m = size(a);
  x = x / norm(x);
  
  if sigma == inf
    sigma = x' * a * x;
  end

  res = norm((a - sigma*speye(m))*x);
  
  eigval_iterates = [sigma];  % save the approximations for debugging,
  eigvec_iterates = [x];      % plotting or convergence analysis
  residuals = [res];
  
  iterations = 0;

  while res >= tolerance
    x = (a - sigma*speye(m) + res*1i*(speye(m) - x*x'))\x;
    x = x / norm(x);
    sigma = (x' * a * x);
    res = norm((a - sigma*speye(m))*x);
    
    eigval_iterates = [eigval_iterates, sigma]; %append current iterates
    eigvec_iterates = [eigvec_iterates, x];     % to list of approx.
    resold = residuals(end);
    residuals = [residuals, res];
    
    iterations = iterations + 1;
  end
  x = x / norm(x);
  sigma = real(sigma);
end
