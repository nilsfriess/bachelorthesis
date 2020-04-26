dimension = 10;
angle = pi/2;

% Define pdf for phi ~ sin^(n-2)(phi)
samples = 1000;
xs = linspace(0, deg2rad(angleDeg), samples);
pdf = sin(xs) .^ (N-2);

% normalise pdf so that it integrates to 1
pdf = pdf / sum(pdf);

% construct cdf
cdf = cumsum(pdf);
[cdf, mask] = unique(cdf);
xs = xq(mask);

% Perform actual inversion sampling
phi = interp1(cdf, xq, rand());

% Convert from spherical -> cartesian
x1 = cos(phi);

% Sample remaining coordinates from N-1-dimensional sphere
x2_n = randn(N-1, 1);
x2_n = x2_n / norm(x2_n);
x = [x1;x2_n];
x = x / norm(x);