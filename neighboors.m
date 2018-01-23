function nei = neighboors(graph)
    %compute the adjacency list of the graph
    [n,~] = size(graph.Nodes);
    nei = cell([1,n]);
    g1 = graph.Edges.EndNodes(:,1);
    g2 = graph.Edges.EndNodes(:,2);
    parfor i =1:n
        a = g1(g2==i);
        b = g2(g1==i);
        neig = transpose(cat(1,a,b));
        nei{i} = neig;
    end
    end