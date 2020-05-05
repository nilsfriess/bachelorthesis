function z = laplace_eigfun_exact(k,l,x,y)
    % Taken from http://people.inf.ethz.ch/arbenz/ewp/Slides/slides1.pdf
    z = sin(k*pi.*x).*sin(((2*l -1)/2)*pi.*y);
end

