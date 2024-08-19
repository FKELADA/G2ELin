function [paths, edges, weights] = allpaths(G, startnode, endnode, visited, path, edge_path, weight_path)
    if nargin < 3
        error('Not enough input arguments.');
    end

    if nargin == 3
        visited = [];
        path = [];
        edge_path = [];
        weight_path = [];
    end

    visited = [visited, startnode];
    path = [path, startnode];

    if startnode == endnode
        paths = {path};
        edges = {edge_path};
        weights = {weight_path};
    else
        paths = {};
        edges = {};
        weights = {};

        neighbors = G.neighbors(startnode);

        for i = 1:length(neighbors)
            neighbor = neighbors(i);
            edge = findedge(G, startnode, neighbor);
            weight = G.Edges.Weight(findedge(G, startnode, neighbor)); % assuming weights are stored in the 'Weight' property

            if ~ismember(neighbor, visited)
                newpaths = allpaths(G, neighbor, endnode, visited, path, [edge_path, edge], [weight_path, weight]);
                paths = [paths, newpaths];
                
                % Find the edges and weights between consecutive nodes in the path
                new_edges = cellfun(@(p) find_edges_between_nodes(G, p), newpaths, 'UniformOutput', false);
                new_weights = cellfun(@(p) find_weights_between_nodes(G, p), newpaths, 'UniformOutput', false);
                edges = [edges, new_edges];
                weights = [weights, new_weights];
            end
        end
    end
end

function edge_path = find_edges_between_nodes(G, path)
    edge_path = [];
    for i = 1:length(path)-1
        edge = findedge(G, path(i), path(i+1));
        edge_path = [edge_path, edge];
    end
end

function weight_path = find_weights_between_nodes(G, path)
    weight_path = [];
    for i = 1:length(path)-1
        weight = G.Edges.Weight(findedge(G, path(i), path(i+1))); % assuming weights are stored in the 'Weight' property
        weight_path = [weight_path, weight];
    end
end
