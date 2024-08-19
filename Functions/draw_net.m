function [G] = draw_net(Y_network, Y_line)
s = Y_line(:,2);
t = Y_line(:,3);
R = round(Y_line(:, 4), 4);
X = round(Y_line(:, 5), 4);
Line_RX = R ./ X;
Line_z = sqrt(R.^2+X.^2);
weights = round(Line_z,3);

% Define colors based on node type
none_color = [0 0 0];
inf_bus_color = [0.64 0.08 0.18];
grid_forming_color = [0 0.8 0];
grid_following_color = [0 0 0.8];
sync_machine_color = [0.85 0.33 0.1];

% Create a vector of node colors and shapes based on the Y_network matrix
node_colors = zeros(size(Y_network, 1), 3);
node_shapes = repmat({'o'}, size(Y_network, 1), 1);

for i = 1:size(Y_network, 1)
    node_type = Y_network(i, 2);
    if node_type == 0
        node_colors(i, :) = none_color;
    elseif node_type == 1
        node_colors(i, :) = inf_bus_color;
    elseif node_type == 2
        node_colors(i, :) = grid_forming_color;
    elseif node_type == 3
        node_colors(i, :) = grid_following_color;
    elseif node_type == 4
        node_colors(i, :) = sync_machine_color;
    end
    
    load_status = Y_network(i, 3);
    if load_status == 1
        node_shapes{i} = 's';
    end
end

% Create the graph object and plot it with customized node colors and shapes
G = graph(s, t, weights);
% ind = 1:2:size(Y_network, 1);
% Gnew = gsp_kron_reduction(G,ind);
figure
plot(G,'EdgeLabel', G.Edges.Weight, 'NodeColor', node_colors,'Marker',node_shapes,'MarkerSize', 7,'LineWidth',2);

end