function [string_values] = modal_analysis_sens(plot_order,A_tot0,B_tot0,C_tot0,D_tot0,stateNames,inputNames,outputNames,mode_sensitivity,mode_participation,mode_mshape,n_exc_ln_ld,start_from_eig,end_eig,StepResp_Data,FreeResp_Data,init_V,xz,yz,Balance_option)
sens_plot = plot_order(1);
part_matrix_plot = plot_order(2);
single_mode_part = plot_order(3);
eig_plot = plot_order(4);
mode_shape = plot_order(5);
free_motion = plot_order(6);
step_resp = plot_order(7);
%% Modal Analysis:
% 1) Eigenvalue calc:
[Vu,Du,Wu] = eig(A_tot0,Balance_option); % obtain eigenvalues and eigenvectors
% columns of Vu are the UNSORTED right eigenvectors
% columns of Wu are the UNSORTED left eigenvectors
% diagonal of Du are the UNSORTED elements of the eigenvalues
% Sorted order:
[d,ind] = sort(diag(Du),'descend');
D = Du(ind,ind);
eigenVal = d; %eigenVal are the sorted eigenvalues 
% V = Vu(:,ind);
% W = Wu(:,ind);
V = Vu(:,ind);
W = Wu(:,ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Eigenvectors:
% [V,D,W] = eig(A_tot0,Balance_option); % obtain eigenvalues and eigenvectors
% Eigenvectors
% columns of V are the right eigenvectors
% columns of W are the left eigenvectors
coef = diag(W'*V); % for normalization of eigenvectors
rightV = V;
leftV = zeros(size(A_tot0));
for i = 1:length(A_tot0)
    leftV(i,:) = 1/coef(i)*(W(:,i)'); % normalization of left eigenvectors
end
% columns of rightV are the right eigenvectors
% rows of leftV are the left eigenvectors
% rightV and leftV are orthonormalized -> leftV*rightV = I
% ATTENTION : leftV is not like Matlab W, it is actually W' normalized;
% the ~rows~ of leftV are eigenvectors, not the columns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Sensitivity Matrix:
% What is the mode you want to study its sensitivity?
if sens_plot == 1
mode = mode_sensitivity; 
sensMatrix = zeros(length(A_tot0),length(A_tot0),length(A_tot0));
for i = 1:length(A_tot0)
    for k = 1:length(A_tot0)
        for j = 1:length(A_tot0)
            sensMatrix(k,j,i) = (leftV(i,k)*rightV(j,i))/(leftV(i,:)*rightV(:,i));
        end
    end
end

Anum_ones = A_tot0./A_tot0;
Anum_ones(isnan(Anum_ones))=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4) Participation Matrix:
% Heat map of the participation factors of all states in all modes:
for i = 1:numel(eigenVal)
    xlabels{i} = sprintf('\\lambda_{%d}',i);
end
ylabels = stateNames;
partMatrix = zeros(size(A_tot0));
den = zeros(length(A_tot0));
for j = 1:length(A_tot0)
    sum = 0;
    for k = 1:length(A_tot0)
        sum = sum + abs(leftV(j,k))*abs(rightV(k,j));
    end
    den(j) = sum;
end
for i = 1:length(A_tot0)
    for j = 1:length(A_tot0)
        partMatrix(i,j) = abs(rightV(i,j))*abs(leftV(j,i))/den(j);
        % rightV(i,j) is the activity of ith state in jth mode
        % leftV(j,i) is the contribution of this activity to jth mode
        % P(i,j) is the relative participation of ith state in jth mode
        % each row of P corresponds to a state
        % each column of P corresponds to a mode
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) EigenValues Analysis & plot:
nTopStates = 3;
% Initialize variables for storing results
num_eigs = length(eigenVal);
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
    damp_ratio(i) = -real(eigenVal(i))/abs(eigenVal(i));
    undamped_freqs(i) = abs(eigenVal(i))/(2*pi); %in Hz
    damped_freqs(i) = undamped_freqs(i) * sqrt(1 - damp_ratio(i)^2);
end


% Plotting:
zita_min    = 0.05; %minimum damping ratio
zita_lim    = min(real(eigenVal(:)));   %trace the minimum damping axe with the slope of the minimum damping between these frequency range




end