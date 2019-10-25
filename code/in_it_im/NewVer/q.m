function pot = q(mid)

xy = 40*mid-20;
r2 = xy(:,1).^2 + xy(:,2).^2;

pot = - 5*exp(-r2) + cos(xy(:,1)) + cos(xy(:,2));

end