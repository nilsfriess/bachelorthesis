function x = random_in_cone(direction, angleDeg)
    direction = direction / norm(direction);
    [N,~] = size(direction);
    
    % Define pdf
    samples = 5000;
    xs = linspace(0, deg2rad(angleDeg), samples);
    pdf = sin(xs) .^ (N-2);
    
    % normalise pdf so that it integrates to 1
    pdf = pdf / sum(pdf);
    
    % build cdf
    cdf = cumsum(pdf);
    xq = xs;
    [cdf, mask] = unique(cdf);
    xq = xq(mask);
    
    % perform actual inversion sampling
    phi = interp1(cdf, xq, rand());
    
    % spherical -> cartesian
    x1 = cos(phi);
    
    % sample remaining coordinates from N-1-dimensional sphere
    x2_n = randn(N-1,1);
    x2_n = x2_n / norm(x2_n);
    x = [x1;x2_n];
    x = x / norm(x);
end