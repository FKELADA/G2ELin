function [string_values] = sumTable(eigA,nTopStates,partMatrix,ylabels,xlabels)

    % Initialize variables for storing results
    num_eigs = length(eigA);
    top_states = nan(num_eigs, nTopStates); % top "nTopStates" associated states for each eigenvalue
    damped_freqs = zeros(num_eigs, 1); % damped frequency for each eigenvalue
    undamped_freqs = zeros(num_eigs, 1); % undamped frequency for each eigenvalue
    dampings = zeros(num_eigs, 1); % damping for each eigenvalue pair

    % Loop over sorted eigenvalues
    for i = 1:num_eigs
        % Find top "nTopStates" associated states for each eigenvalue
        input_matrix = abs(partMatrix);
        vector_matrix = ylabels;
        [num_rows, num_cols] = size(input_matrix);
        highest_values = zeros(nTopStates, 2);
        % Loop through each column of the input matrix
        % Get the values and row indices for the current column
        [col_values, row_indices] = sort(input_matrix(:, i), 'descend');
        % Store the highest 3 values and their row indices in the pre-allocated matrix
        for k = 1:nTopStates
            highest_values(k, :) = [col_values(k), row_indices(k)];
        end
        
        % Get the linear indices of the highest values in the vector matrix
        linear_indices = sub2ind([num_rows, 1], highest_values(:, 2));
        % Get the strings at the cell indices of the highest values in the vector matrix
        string_values(i,:) = vector_matrix(linear_indices);

        % Calculate damped frequency and undamped frequency for each eigenvalue
        damp_ratio(i) = -real(eigA(i))/abs(eigA(i));
        undamped_freqs(i) = abs(eigA(i))/(2*pi); %in Hz
        damped_freqs(i) = undamped_freqs(i) * sqrt(1 - damp_ratio(i)^2);
    end



    % Display results in a table
    T = table(xlabels', real(eigA),imag(eigA), string_values(:,1), string_values(:,2), string_values(:,3), undamped_freqs, damped_freqs, damp_ratio'*100);
    T.Properties.VariableNames = {'Eigenvalue No.','\Re(\lambda)','\Im(\lambda)', 'State1', 'State2' , 'State3', 'UndampedFrequency (Hz)', 'DampedFrequency (Hz)', 'DampingPercentage (%)'};
    disp(T);

    data = [xlabels', num2cell(real(eigA)),num2cell(imag(eigA)), string_values(:,1), string_values(:,2), string_values(:,3), num2cell(undamped_freqs), num2cell(damped_freqs), num2cell(damp_ratio'*100)] ;
    CN   = {'Eigenvalue No.','$\Re(\lambda)$','$\Im(\lambda)$', 'State1', 'State2', 'State3', 'UndampedFrequency (Hz)', 'DampedFrequency (Hz)', 'DampingPercentage (%)'};
    t = uitable(uifigure,'Data', data, 'ColumnName', CN);

end