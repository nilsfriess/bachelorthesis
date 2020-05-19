function [x, sigma, iterations, eigval_iterates, eigvec_iterates, residuals] = classic_rqi(a, x, sigma, tolerance)
% CLASSIC_RQI   Computes an eigenpair of a using the Rayleigh quotient
% iteration
%   [x, sigma, iterations, eigval_iterates, eigvec_iterates] = CLASSIC_RQI(a, x, sigma) 
%   Uses sigma as an initial guess for the eigenvalue and x as an initial 
%   guess for the eigenvector of the matrix a. If, however, sigma == inf, 
%   the Rayleigh quotient of x w.r.t. a is used for the initial eigval 
%   guess. The final approximated eigenvector and -value are stored in x
%   and sigma, repsectively. The results after every iteration are stored 
%   in the vectors eigval_iterates and eigvec_iterates. The iteration is 
%   stopped if the residual is less than 1e-9.
%
%   [x, sigma, iterations, eigval_iterates, eigvec_iterates] = CLASSIC_RQI(a, x, sigma, tol) 
%   Same as above except that the computation is stopped when the residual 
%   is less than tol.

    if nargin < 3
        error("Too few arguments");
    elseif nargin < 4
        % if no tolerance provided, use 10e-9
        tolerance = 1e-11; 
    end

    m = size(a);
    x = x / norm(x); % normalise the eigenvector guess
    
    if sigma == inf
       sigma =  x' * a * x;
    end
    
    eigval_iterates = [sigma];  % save the approximations for debugging,
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