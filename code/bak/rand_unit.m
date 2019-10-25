function v = rand_unit(m)
    v = random("unif",-1,1,[m,1]);
    v = v / norm(v);
end