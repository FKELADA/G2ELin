function [results_matrix] = SCR(G_graph,Y_network,Y_DER,Sb_MW,Un_up,Un_line,Zb_line)

Xd = Y_DER(1,13);
% Initialize matrix to store results
results_matrix = zeros(size(Y_network, 1), 2); % Two columns: Node and SCR Value
results_Isc = zeros(size(Y_network, 1), 2);
% Display header for the table
disp('Node   SCR Value');
n_nodes = size(Y_network, 1);
n_DER = size(Y_DER, 1);
% Loop over all nodes
for x = 1:n_nodes
    % Initialize equivalent admittance for each iteration
    M_1Path = Sb_MW/Xd;
    M_AllPaths = M_1Path; 

    % Find all possible paths from node 1 to node x
    [paths,edges,weights] = allpaths(G_graph, n_nodes-n_DER+1, x);

    % Calculate equivalent admittance (sum of admittances along each path)

    for i = 1:size(paths{1},1)
        path = paths{i};
        edge = edges{i};
        weight = weights{i};        
        for j = 1:size(edges{1},2)
            M_1Edge = (Un_line)^2./(Zb_line*weight(j));
            M_1Path = M_1Path*M_1Edge/(M_1Path+M_1Edge);
        end
        M_AllPaths = M_1Path;
    end


    % Calculate Short Circuit Current (Isc)
    if x == n_nodes-n_DER+1
        V_S = Un_up;
    else
        V_S = Un_line;
    end
    Isc = M_AllPaths / (sqrt(3) * V_S);
    In = Sb_MW / (sqrt(3) * Un_up);
    SCR_values = Isc / In;
    
    % Display node index and SCR value
    fprintf('%4d   %.4f\n', x, SCR_values);

    % Store node index and SCR value in the results matrix
    results_matrix(x, :) = [x, SCR_values];
    results_Isc(x,:) = [x,Isc];
end



end