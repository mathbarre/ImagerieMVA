% function nei=neighboors(i,graph)
% a = graph.Edges.EndNodes(graph.Edges.EndNodes(:,2)==i,1);
% b = graph.Edges.EndNodes(graph.Edges.EndNodes(:,1)==i,2);
% nei = transpose(cat(1,a,b));
% end

function nei = neighboors(graph)
[n,~] = size(graph.Nodes);
nei = cell([1,n]);
for i =1:n
    a = graph.Edges.EndNodes(graph.Edges.EndNodes(:,2)==i,1);
    b = graph.Edges.EndNodes(graph.Edges.EndNodes(:,1)==i,2);
    neig = transpose(cat(1,a,b));
    nei{i} = neig;
end
end