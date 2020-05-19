function [x, sigma, iterations, eigval_iterates, eigvec_iterates, residuals] = complex_rqi2combined(a,x,sigma,gamma,tolerance)
% complex_rqi   Computes an eigenpair of a using the complex Rayleigh 
% quotient iteration
%
%   See CLASSIC_RQI for description of return values. The parameters a, x, 
%   sigma and tol are also described there.

  if nargin < 4
    error("Too few arguments");
  elseif nargin < 5
    tolerance = 1e-11;
  end
  m = size(a);
  x = x / norm(x);
    
  if sigma == inf
    sigma = x' * a * x;
  end

  res = norm((a - sigma*speye(m))*x);
  
  if gamma == inf
      if res > 1
          gamma = res;
      else
          gamma = res^2;
      end
  end
  
  eigval_iterates = [sigma];  % save the approximations for debugging,
  eigvec_iterates = [x];      % plotting or convergence analysis
  residuals = [res];
  
  iterations = 0;
    
  while res >= tolerance
    x = (a - (sigma + gamma*1i)*speye(m))\x;
    x = x / norm(x);
    sigma = x'*a*x;
        
    res = norm((a - sigma*speye(m))*x);
    
    if res >= 1
       gamma = res;
    else
       gamma = res^2;
    end
    eigval_iterates = [eigval_iterates, sigma]; %append current iterates
    eigvec_iterates = [eigvec_iterates, x];     % to list of approx.
    residuals = [residuals, res];
        
    iterations = iterations + 1;
  end
  
  x = real(x) / norm(real(x));
  sigma = x'*a*x;
  residuals(end) = norm((a - sigma*speye(m))*x); 
end
