% Function to set up the initial grid.
%
% You can replace the mesh by whatever you like (satisfying
% the rules set out at the beginning of Section 2.2). Just
% edit the function if you want to change the mesh. Use 
% the function "visualise" to check that you have the mesh 
% which you want.

function [nodes,elements] = initial_mesh(topright)

nodes = [[0.0 0.0 1];...
         [1.0 0.0 1];...
         [1.0 topright 1];...
         [0.0 1.0 1]];    % The entry in the third column specifies 
                          % whether the node is a boundary node or
                          % an interior one (1=boundary, 0=interior)

elements = [[1 2 3 ];...
            [1 3 4 ]];
