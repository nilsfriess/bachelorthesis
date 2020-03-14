function [x, sigma, iterations, eigval_iterates, eigvec_iterates, residuals] = complex_rqi_AMG(a,x,sigma,gamma,tolerance)
% complex_rqi   Computes an eigenpair of a using the complex Rayleigh 
% quotient iteration
%
%   See CLASSIC_RQI for description of return values. The parameters a, x, 
%   sigma and tol are also described there. gamma is the initial imaginary
%   shift. If gamma == Inf, the residual is used as the initial i-shift.

    if nargin < 4
        error("Too few arguments");
    elseif nargin < 5
        % if no tolerance provided, use 10e-9
        tolerance = 10e-9;
    end
    m = size(a);
    x = x / norm(x);
    
    I = speye(m);
    
    if sigma == inf
       sigma =  x' * a * x;
    end
    
    if gamma == inf
       gamma = norm((a - sigma*I)*x);
    end
    res = gamma;
    
    eigval_iterates = [sigma];  % save the approximations for debugging,
    eigvec_iterates = [x];      % plotting or convergence analysis
    residuals = [res];

    A = a - (sigma + res*res*1i)*I;
    
    iterations = 0;
    
    options.isdefinite = 1;
    options = AMGinit(A, options);
    [PREC, options] = AMGfactor(A, options)
        
    while res >= tolerance
        A = a - (sigma + res*res*1i)*I;
        options.restol = 10e-5;
        [x,options] = AMGsolver(A, PREC, options, x, x); % last is initial guess

        x = x / norm(x);
        % x = real(x);
        sigma = x' * a * x;
        res = norm((a - sigma*I)*x);
        
        eigval_iterates = [eigval_iterates, sigma]; %append current iterates
        eigvec_iterates = [eigvec_iterates, x];     % to list of approx.
        residuals = [residuals, res];
        
        iterations = iterations + 1;
    end
    x = real(x);
    x = x / norm(x);
    sigma = real(sigma);
end