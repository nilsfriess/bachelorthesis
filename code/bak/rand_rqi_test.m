function [reseigs, resinitv, diffv] = rand_rqi_test(m, n, shift, exactv)
    a = generate_test_matrix(m+1);
    reseigs = zeros(n,1);
    resinitv = zeros(n,m);
    exactv = exactv / norm(exactv);
    diffv = zeros(n,1);
       
    for i = 1:n
        v = rand_unit(m);
        [~, ee] = classic_rqi(a, v, shift);
        reseigs(i) = ee;
        resinitv(i,:) = v;
        diffv(i) = v'*exactv;
    end
    diffv = rad2deg(acos(diffv));
end