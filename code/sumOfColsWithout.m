function s = sumOfColsWithout(a,col)
    [r,~] = size(a);
    s = zeros(r,1);
    
    for i = 1:r
        % sum all columns expect for col
        if r ~= col
            s = s + a(:,i);
        end
    end
end

