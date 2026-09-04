clc; clear;
% Sample Y_network matrix (Node admittance matrix)
Y_network = [
    1   0   0   0; % Node 1: Slack node
    0   1   0   0; % Node 2: PV node
    0   0   1   0; % Node 3: PQ node
];

% Sample Y_line matrix (Line admittance matrix)
Y_line = [
    1   -0.5    -0.5; % Line 1 between Node 1 and Node 2
    -0.5    1   -0.5; % Line 2 between Node 2 and Node 3
    -0.5   -0.5   1;  % Line 3 between Node 1 and Node 3
];

% Number of nodes
num_nodes = size(Y_network, 1);

% Initialize SCR array
SCR_values = zeros(num_nodes, 1);

% Iterate over each node
for i = 1:num_nodes
    % Extract Y_bus for the specific node
    Y_bus = Y_network(i, i);

    % Extract the lines connected to the node
    connected_lines = find(Y_line(:, i));

    % Sum the line admittances
    for j = 1:length(connected_lines)
        line_index = connected_lines(j);
        Y_bus = Y_bus + Y_line(line_index, line_index);
    end

    % Calculate SCR for the node
    SCR_values(i) = abs(Y_bus) / real(Y_bus);
end

% Display SCR values for each node
disp('Short Circuit Ratios (SCR) for each node:');
disp(SCR_values);
