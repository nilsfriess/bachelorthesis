function node = visited(edges,node1,node2)

[nedges,m]  = size(edges);

node = 0;
for i=1:nedges
    
    if ((edges(i,1) == node1 & edges(i,2) == node2) |...
        (edges(i,1) == node2 & edges(i,2) == node1))
        node = i;
        return
    end

end

