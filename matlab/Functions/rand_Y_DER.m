function [Y_DER] = rand_Y_DER(Y_DER_i)

% Step 1: Generate random values for the second column
% second_column = zeros(size(Y_DER_i, 1),1);
% 
% % Define probabilities
% val_2_4 = 0;
% 
% % Generate random values
% for i = 1:numel(second_column)
%     rand_val = randi([2 4],1);
%     if rand_val == 2 || rand_val == 4
%         val_2_4 = val_2_4 + 1;
%         perc_2_4 = val_2_4/14;
%         if perc_2_4 <= 0.5
%             second_column(i) = rand_val;
%         else 
%             second_column(i) = 3;
%         end
%     else 
%         second_column(i) = 3;
%     end
% end
second_column = randi([2, 4], size(Y_DER_i, 1), 1);

% Ensure at least one value in the second column is 2 or 4
while ~any(second_column == 2) && ~any(second_column == 4)
    second_column = randi([2, 4], size(Y_DER_i, 1), 1);
end
Y_DER_i (:,2) = second_column;

% Step 2: Generate random values for the fourth column
for i_n = 1:14
    if Y_DER_i(i_n,2) == 2 %it is a GFM
        Y_DER_i(i_n,4) = 2;
    elseif Y_DER_i(i_n,2) == 3 %it is a GFL
        Y_DER_i(i_n,4) = 3;
    elseif Y_DER_i(i_n,2) == 4 %it is a SM
        Y_DER_i(i_n,4) = 2;
    end
end

% Select a row where the second column is 2 or 4 and set its fourth column to 1
indices_2_or_4 = find(second_column == 2 | second_column == 4); % Indices where second column is 2 or 4
if ~isempty(indices_2_or_4)
    % Randomly select a row index from the indices where fourth column is 2
    random_index = indices_2_or_4(randi(length(indices_2_or_4)));
    % Set the selected row in column 4 to 1
    Y_DER_i(random_index,4) = 1;
end 


% Step 3: Generate random values for the seventh column with the sum constraint
max_sum = 3; % Constraint on the total sum
seventh_column = zeros(size(Y_DER_i, 1), 1);
rand_values = rand(size(Y_DER_i, 1), 1) * 0.4 + 0.1; % Generate random values between 0.1 and 0.5

% Calculate the scaling factor to ensure the sum constraint
scaling_factor = max_sum / sum(rand_values);

% Apply the scaling factor to the random values
seventh_column = rand_values * scaling_factor;

% Update the seventh column
Y_DER_i(:, 7) = seventh_column;

second_column = Y_DER_i(:,2);
% Find the indices of the rows where the second column is 2, 4, and 3
index_2 = find(second_column == 2);

Y_DER = [];
% Move the rows where the second column is 2 to the top of the matrix
if ~isempty(index_2)
    temp_row_2 = Y_DER_i(index_2, :);
    Y_DER = temp_row_2;
end

index_4 = find(second_column == 4);
% Move the rows where the second column is 4 to the top of the matrix
if ~isempty(index_4)
    temp_row_4 = Y_DER_i(index_4, :);
    Y_DER = [Y_DER; temp_row_4];
end

index_3 = find(second_column == 3);
if ~isempty(index_3)
    temp_row_3 = Y_DER_i(index_3, :);
    Y_DER = [Y_DER; temp_row_3];
end
second_column = Y_DER(:,2);
fourth_column = Y_DER(:,4);
% Find the index of the row where the fourth column is 1
index_1 = find(fourth_column == 1);

% Move the rows where the fourth column is 1 to the top of the matrix
if ~isempty(index_1)
    % Remove the row where fourth column is 1 from the original matrix
    temp_row = Y_DER(index_1, :);
    Y_DER(index_1, :) = [];
    % Concatenate the removed row at the beginning of the matrix
    Y_DER = [temp_row; Y_DER];
end
