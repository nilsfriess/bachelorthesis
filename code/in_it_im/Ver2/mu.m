function a = mu(X)

% Takes three nodes of a triangle, listed as the rows 
% of a matrix X and finds the area of this triangle

a =  abs(det([ones(3,1),X])/2);
