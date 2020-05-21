function [x, sigma, iterations, eigval_iterates, eigvec_iterates, residuals] = classic_rqi(a, x, sigma, tolerance)
% CLASSIC_RQI   Computes an eigenpair of the matrix 'a' using the 
% Rayleigh quotient iteration
%   [x, sigma, iterations, eigval_iterates, eigvec_iterates] = CLASSIC_RQI(a, x, sigma) 
%   Uses sigma as an initial guess for the eigenvalue and x as an initial 
%   guess for the eigenvector of the matrix a. If sigma == inf, 
%   the Rayleigh quotient of x w.r.t. a is used for the initial eigval 
%   guess. The final approximated eigenvector and eigenvalue are stored in 
%   x and sigma, respectively. The remaining return arguments contain the
%   number of iterations and the sequences of approximations of the 
%   eigenvalue, eigenvectors and residuals, respectively.
%
%   [x, sigma, iterations, eigval_iterates, eigvec_iterates] = CLASSIC_RQI(a, x, sigma, tol) 
%   Same as above except that the computation is stopped when the residual 
%   is less than tol.

    if nargin < 3
        error("Too few arguments");
    elseif nargin < 4
        % if no tolerance provided, use 10^(-9)
        tolerance = 1e-9; 
    end

    m = size(a);
    x = x / norm(x);
    
    if sigma == inf
       sigma =  x' * a * x; % Compute Rayleigh Quotient
    end
    
    eigval_iterates = [sigma];  % Save the approximations for debugging,
    eigvec_iterates = [x];      % plotting or convergence analysis
    
    res = norm((a - sigma*speye(m))*x);
    residuals = [res];
    
    iterations = 0;
    while res >= tolerance
       x = (a - sigma * speye(m)) \ x;  % solve linear system
       x = x / norm(x);                 % normalise iterate
       sigma = x' * a * x;              % compute Rayleigh quotient
       
       res = norm((a - sigma*speye(m))*x);
       
       eigval_iterates = [eigval_iterates, sigma]; %append current iterates
       eigvec_iterates = [eigvec_iterates, x];     % to list of approx.
       residuals = [residuals, res];       
       
       iterations = iterations + 1;
    end
end