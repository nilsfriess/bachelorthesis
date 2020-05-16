function qs = convergence_rate(vals, limit)
    if length(vals) <= 3
        error("Too few values");
    end
%     qs = zeros(length(vals)-1,1);
%     for i = 2 : length(vals)
%        a = 1 / log(2) * log( norm( (vals(i-1) - limit) / ( vals(i) - limit) ));
%        qs(i-1) = a;
%     end

    % From: https://en.wikipedia.org/wiki/Rate_of_convergence
    qs = zeros(length(vals)-2,1);
    for i = 3:length(vals)-1
        j = i-2;
        qs(j) = log(norm((vals(i+1) - vals(i))/(vals(i)-vals(i-1))));
        qs(j) = qs(j) / log(norm((vals(i) - vals(i-1))/(vals(i-1) - vals(i-2))));
    end
end
