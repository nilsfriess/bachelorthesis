function convergence_rate(vals, limit)
    if length(vals) <= 3
        error("Too few values");
    end
    for i = 2 : length(vals)
       a = 1 / log(2) * log( norm( (vals(i-1) - limit) / ( vals(i) - limit) ));
       disp(a);
    end
end
