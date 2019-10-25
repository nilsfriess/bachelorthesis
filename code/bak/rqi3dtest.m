function rqi3dtest()
    tiledlayout(1,2)

    nexttile
    [reseigs,ri,~] = rand_rqi_test_complex(3, 200, 30, 10, [-0.7071;0;0.7071]);
    [m,~] = size(ri);
    for i = 1:m
       vec = ri(i,:);
       reseigs(i)
       if abs(reseigs(i) - 32) < 10e-2
           color = "blue";
       else
           color = "red";
           endQ
       quiver3(0,0,0,vec(1), vec(2), vec(3), 'color', color);
        
       hold on
    end
    
    nexttile
    [reseigs,ri,~] = rand_rqi_test(3, 200, 30, [-0.7071;0;0.7071]);
    [m,~] = size(ri);
    for i = 1:m
       vec = ri(i,:);
       reseigs(i)
       if abs(reseigs(i) - 32) < 10e-2
           color = "blue";
       else
           color = "red";
       end
       quiver3(0,0,0,vec(1), vec(2), vec(3), 'color', color);
        
       hold on
    end
end

