function [] = modal_analysis(plot_order,A_tot0,B_tot0,C_tot0,D_tot0,stateNames,inputNames,outputNames,mode_sensitivity,mode_participation,mode_mshape,n_exc_ln_ld,start_from_eig,end_eig,StepResp_Data,FreeResp_Data,init_V,xz,yz,Balance_option,latexTab,selected_modes)
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
max_sens = cell(num_eigs, 8);

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
    pmode = sort(partMatrix(:,i),'descend');
    parmode(i,1:3) = pmode(1:3,1);
    
    if sens_plot == 1
    % max sens matrix A indices
    % Get the matrix slice at z
    matrix_slice = abs(sensMatrix(:,:,i)).*Anum_ones;
    % Sort the matrix slice in descending order
    [sorted_values, sorted_indices] = sort(matrix_slice(:), 'descend');
    % Get the top 5 values and their indices
    top_8_values = sorted_values(1:8);
    top_8_indices = sorted_indices(1:8);
    % Convert indices to subscripts (x,y)
    [x_ind, y_ind] = ind2sub(size(matrix_slice), top_8_indices);
    % Format indices as strings and store them in max_sens array
    for j = 1:8
        max_sens{i, j} = sprintf('A_{%d_{%d}}(%0.3f)', x_ind(j), y_ind(j),top_8_values(j));
    end
    end
end

if sens_plot == 1
% Display results in a table
T = table(xlabels', real(eigenVal), imag(eigenVal), string_values(:,1), parmode(:,1).*100, string_values(:,2), parmode(:,2).*100, string_values(:,3), parmode(:,3).*100, undamped_freqs, damped_freqs, damp_ratio'*100, max_sens);
T.Properties.VariableNames = {'Eigenvalue No.','\Re(\lambda)','\Im(\lambda)', 'State1', 'Part1. (%)', 'State2', 'Part2. (%)', 'State3', 'Part3. (%)', 'UndampedFrequency (Hz)', 'DampedFrequency (Hz)', 'DampingPercentage (%)','A matrix sens. elem.'};
disp(T);

data = [xlabels', num2cell(real(eigenVal)),num2cell(imag(eigenVal)), string_values(:,1), parmode(:,1).*100, string_values(:,2), parmode(:,2).*100, string_values(:,3), parmode(:,3).*100, num2cell(undamped_freqs), num2cell(damped_freqs), num2cell(damp_ratio'*100),max_sens] ;
CN   = {'Eigenvalue No.','$\Re(\lambda)$','$\Im(\lambda)$', 'State1', 'Part1. (%)', 'State2', 'Part2. (%)', 'State3', 'Part3. (%)', 'UndampedFrequency (Hz)', 'DampedFrequency (Hz)', 'DampingPercentage (%)','A matrix sens. elem.'};
t = uitable(uifigure,'Data', data, 'ColumnName', CN);
else
    % Display results in a table
T = table(xlabels', real(eigenVal), imag(eigenVal), string_values(:,1), parmode(:,1).*100, string_values(:,2), parmode(:,2).*100, string_values(:,3), parmode(:,3).*100, undamped_freqs, damped_freqs, damp_ratio'*100);
T.Properties.VariableNames = {'Eigenvalue No.','\Re(\lambda)','\Im(\lambda)', 'State1', 'Part1. (%)', 'State2', 'Part2. (%)', 'State3', 'Part3. (%)', 'UndampedFrequency (Hz)', 'DampedFrequency (Hz)', 'DampingPercentage (%)'};
disp(T);

data = [xlabels', num2cell(real(eigenVal)),num2cell(imag(eigenVal)), string_values(:,1), parmode(:,1).*100, string_values(:,2), parmode(:,2).*100, string_values(:,3), parmode(:,3).*100, num2cell(undamped_freqs), num2cell(damped_freqs), num2cell(damp_ratio'*100)] ;
CN   = {'Eigenvalue No.','$\Re(\lambda)$','$\Im(\lambda)$', 'State1', 'Part1. (%)', 'State2', 'Part2. (%)', 'State3', 'Part3. (%)', 'UndampedFrequency (Hz)', 'DampedFrequency (Hz)', 'DampingPercentage (%)'};
t = uitable(uifigure,'Data', data, 'ColumnName', CN);
end
% Plotting:
zita_min    = 0.05; %minimum damping ratio
zita_lim    = min(real(eigenVal(:)));   %trace the minimum damping axe with the slope of the minimum damping between these frequency range


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting:

% Plotting the sensitivity matrix:
if sens_plot == 1
    f = figure;
    f.Position = [50 50 600 600];  % Max for IEEE paper (687 x 844 pixels)
    thisBar = bar3(abs(sensMatrix(:,:,mode)).*Anum_ones);
    for k = 1:length(thisBar)
        zdata = thisBar(k).ZData;
        thisBar(k).CData = zdata;
        thisBar(k).FaceColor = 'interp';
    end
    colormap(flipud(bone));
    box on;
    thisC = colorbar();
    thisPosC = get(thisC,'Position');
    thisPosC(1) = .93;
    thisPosC(3) = .013;
    set(thisC, 'Position', thisPosC);
    view(0,90);
    set(gca,'Color',[0.94 0.94 0.94]);
    xticklabels([1:length(A_tot0)]);
    xticks([1:length(A_tot0)]);
    yticks([1:length(A_tot0)]);
    xlabel('column (A)');
    ylabel('row (A)');
    yticklabels([1:length(A_tot0)]);
end

% Plotting the participation matrix:
if part_matrix_plot == 1
    figure
    heatmap(xlabels(start_from_eig:end_eig),ylabels(1:n_exc_ln_ld),partMatrix(1:n_exc_ln_ld,start_from_eig:end_eig).*100)
end

% plotting of the participating states in a certain specific mode
% mode_participation
if single_mode_part == 1
    colBlue = [56, 97, 163]/255;
    f = figure;
    f.Position = [50 50 600 600]; % Max for IEEE paper (687 x 844 pixels)
    box on;
    b = bar(abs(partMatrix(:,mode_participation)),'FaceColor',colBlue);
    xticks([1:length(partMatrix)]);
    xticklabels(ylabels);
    xtickangle(90);
    ylabel('Magnitude of part. factor');
end

% % Plotting eigenvalues
% if eig_plot == 1
%     figure 
%     % Group the complex eigenvalues with the same real part and color them similarly
%     real_parts = unique(real(eigenVal(imag(eigenVal)~=0)),'stable'); %the real parts of the poles that have imaginary components ( numel = number of complex poles) 
%     comp_real  = unique(real(eigenVal(imag(eigenVal)==0)),'stable'); %the real parts of the poles that don't have imaginary components ( numel = number of real poles)
%     colors = jet(numel(real_parts)+numel(comp_real));
%     hold on;
%     for i = 1:numel(real_parts)
%         groups = eigenVal(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0);
%         scatter(real(groups), imag(groups)/(2*pi), [], colors(i,:), 'filled');
%         legend_entries{i} = sprintf('[\\lambda_{%d} \\lambda_{%d}] [%s %s %s]', find(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0, 1), find(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0, 1, 'last'),string_values(find(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0, 1),:));
%     end
% 
%     for i = 1:numel(comp_real)
%         gr = eigenVal(real(eigenVal) == comp_real(i) & imag(eigenVal)==0);
%         scatter(real(gr), imag(gr)/(2*pi), [], colors(i+numel(real_parts),:), 'filled');
%         legend_entries{i+numel(real_parts)} = sprintf('\\lambda_{%d} [%s %s %s]', find(real(eigenVal) == comp_real(i),1), string_values(find(real(eigenVal) == comp_real(i),1),:));
%     end
%     hold off;
%     % legend('boxoff')
%     %     xlim([-6500 1]);
%     %     ylim([-400 400]);
%     %     xL = xlim;
%     %     yL = ylim;
%     lgd = legend(legend_entries,'FontSize',6,'AutoUpdate','off','Location','eastoutside');
%     lgd.NumColumns = 2;
%     hold on
%     slope = sqrt((1-(zita_min)^2)/zita_min^2);
%     xz = linspace(zita_lim, 0, 10000);
%     yz = slope * xz/(2*pi);
%     plot(xz, yz, '--r'); % 'b' for blue color
%     hold on
%     plot(xz, -yz, '--r'); % 'b' for blue color
%     hold on
%     slope = sqrt((1-(0.707)^2)/(0.707)^2);
%     xz = linspace(zita_lim, 0, 10000);
%     yz = slope * xz /(2*pi);
%     plot(xz, yz, '--g'); % 'b' for blue color
%     hold on
%     plot(xz, -yz, '--g'); % 'b' for blue color
% 
%     xline(0, 'Color', 'k');  %x-axis
%     yline(0, 'Color', 'k');  %y-axis
%     ylabel('\Im [Hz]');
%     xlabel('\Re');
%     title('Eigenvalues of A');
% end

% Plotting eigenvalues
if eig_plot == 1
    figure 
    % Group the complex eigenvalues with the same real part and color them similarly
    real_parts = unique(real(eigenVal(imag(eigenVal)~=0)),'stable'); %the real parts of the poles that have imaginary components ( numel = number of complex poles) 
    comp_real  = unique(real(eigenVal(imag(eigenVal)==0)),'stable'); %the real parts of the poles that don't have imaginary components ( numel = number of real poles)
    colors = colorcube(numel(real_parts)+numel(comp_real));

    % Set the zoomed-in region
    xlim([xz(1), xz(2)]);
    ylim([yz(1), yz(2)]);

    % Create a logical index for points within the zoomed-in region
    inZoomedRegion = real(eigenVal) >= xz(1) & real(eigenVal) <= xz(2) & imag(eigenVal) >= yz(1)*(2*pi) & imag(eigenVal) <= yz(2)*(2*pi);

    % Create legends only for points within the zoomed-in region
    legend_entries_zoomed = {};
    legend_colors_zoomed = [];
    legend_counter = 1;
    hold on
    for i = 1:numel(real_parts)
        groups = eigenVal(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0 & inZoomedRegion);
        if ~isempty(groups)
            scatter(real(groups), imag(groups)/(2*pi), [], colors(legend_counter,:), 'filled');
            legend_entries_zoomed{legend_counter} = sprintf('[$\\lambda_{%d}$ $\\lambda_{%d}$] [$%s$ $%s$ $%s$]', find(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0 & inZoomedRegion, 1), find(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0 & inZoomedRegion, 1, 'last'), string_values(find(real(eigenVal) == real_parts(i) & imag(eigenVal)~=0 & inZoomedRegion, 1),:));
            legend_colors_zoomed(legend_counter,:) = colors(legend_counter,:);
            legend_counter = legend_counter + 1;
        end
    end

    for i = 1:numel(comp_real)
        gr = eigenVal(real(eigenVal) == comp_real(i) & imag(eigenVal)==0 & inZoomedRegion);
        if ~isempty(gr)
            scatter(real(gr), imag(gr)/(2*pi), [], colors(legend_counter,:), 'filled');
            legend_entries_zoomed{legend_counter} = sprintf('$\\lambda_{%d}$ [$%s$ $%s$ $%s$]', find(real(eigenVal) == comp_real(i) & inZoomedRegion,1), string_values(find(real(eigenVal) == comp_real(i) & inZoomedRegion,1),:));
            legend_colors_zoomed(legend_counter,:) = colors(legend_counter,:);
            legend_counter = legend_counter + 1;
        end
    end
    hold off
    % Add legends only for points within the zoomed-in region
    lgd_zoomed = legend(legend_entries_zoomed, 'FontSize', 10,'AutoUpdate', 'off', 'Location', 'eastoutside', 'Interpreter', 'latex');
    lgd_zoomed.NumColumns = 1;
    hold on
    slope = sqrt((1-(zita_min)^2)/zita_min^2);
    xz = linspace(zita_lim, 0, 10000);
    yz = slope * xz/(2*pi);
    plot(xz, yz, '--r'); % 'b' for blue color
    hold on
    plot(xz, -yz, '--r'); % 'b' for blue color
    hold on
    slope = sqrt((1-(0.707)^2)/(0.707)^2);
    xz = linspace(zita_lim, 0, 10000);
    yz = slope * xz /(2*pi);
    plot(xz, yz, '--g'); % 'b' for blue color
    hold on
    plot(xz, -yz, '--g'); % 'b' for blue color
    
    box on
    xline(0, 'Color', 'k');  %x-axis
    yline(0, 'Color', 'k');  %y-axis
    ylabel('$\Im$ [Hz]', 'Interpreter', 'latex');
    xlabel('$\Re$', 'Interpreter', 'latex');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14);
    title([]);
end

% Plotting the mode shape of a specific mode mode_mshape:
if mode_shape ==1 
    ModeMatrix = zeros(length(A_tot0),2,length(A_tot0));
    for i = 1:length(A_tot0)
        for k = 1:length(A_tot0)
            ModeMatrix(k,1:2,i) = [1 angle(rightV(k,i))];
        end
    end
    nTopStates = 5;
    clear highest_values
    clear string_values
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
        if i == mode_mshape
            index_mshape = linear_indices;
        end
        % Get the strings at the cell indices of the highest values in the vector matrix
        string_values(i,:) = vector_matrix(linear_indices);
    end
    colors = lines(nTopStates);
    figure
    for i = 1:nTopStates
        polarplot([0, ModeMatrix(index_mshape(i),2,mode_mshape)], [0, ModeMatrix(index_mshape(i),1,mode_mshape)],'Color',colors(i,:) , 'LineWidth', 2);
        hold on
    end
    % Add legend with data from the names matrix
    legend(ylabels(index_mshape), 'Location', 'southoutside', 'Orientation', 'horizontal','Position',[0.18214,-0.005714,0.689285,0.059047]);
    title(sprintf('Mode shape of the different states w.r.t mode no. %d',mode_mshape))
end

% Plotting the free motion response of a certain selected state:
if free_motion == 1
    index_state    = find(stateNames == FreeResp_Data(3));
    Tf             = str2double(FreeResp_Data(1));  
    timeStep       = str2double(FreeResp_Data(2));
    C_matrix = leftV * transp(init_V);
    t        = [0:timeStep:Tf];
    for ii = 1:length(eigenVal)
        for i = 1:length(t)
            exp_val(i,ii) = exp(eigenVal(ii)*t(i)); 
        end
    end
    exp_valT = transp(exp_val);
    for i = 1:length(t)
        resp_matrix(:,i) = rightV * (C_matrix.*exp_valT(:,i));
    end
    figure
    title('Free Motion dynamical response of all the states')
    plot(t,resp_matrix(index_state,:))
    xlabel('time')
end


% Plotting the dynamic step response of the system between a certain input-output of the MIMO system:
if step_resp == 1
    step_offset  = str2double(StepResp_Data(1));
    step_Amp     = str2double(StepResp_Data(2));
    Tf           = str2double(StepResp_Data(3));
    stepT        = str2double(StepResp_Data(4));
    input_index = find(inputNames(:)==StepResp_Data(5));
    output_index = find(outputNames(:)==StepResp_Data(6));
    sys               = ss(A_tot0,B_tot0,C_tot0,D_tot0);
    figure
    opt               = stepDataOptions;
    opt.StepAmplitude = step_Amp;
    t                 = 0:stepT:Tf-step_offset;
    y                 = step(sys,t,opt);
    if step_offset == 0
        t     = 0:stepT:Tf;
        y     = step(sys,t,opt);
        plot(t, y(:, output_index, input_index));   %Remember: you have n_outputs and m_inputs such that size(D_tot0) = [n,m] where m -> uk_s (the global inputs)
        title(sprintf('output (%s) after a step change of input (%s) of magnitude %.2f', StepResp_Data(6), StepResp_Data(5), opt.StepAmplitude));
        xlabel('Time');
        ylabel(sprintf('%s', output_index));
    else
        t1   = 0:stepT:(step_offset-stepT);
        t2   = step_offset:stepT:Tf;
        t    = [t1 t2];
        y1   = scaleY*ones(length(t1),1);
        y2   = [y1; scaleY+y(:, output_index, input_index).*scaleY];
        plot(t,y2)
    end

    % sys               = ss(A_tot0,B_tot0,C_tot0,D_tot0);
    % opt               = stepDataOptions;
    % opt.StepAmplitude = 0.1;
    % t                 = 0:1e-5:5;
    % y                 = step(sys,t,opt);
    % save(sprintf('ySMIB_Phasor.mat'),'y');
end

if latexTab == 1
    % Construct LaTeX table content
    latexTable = '\begin{tabular}{|c|c|c|c|c|c|c|c|c|}';
    latexTable = strcat(latexTable, '\hline');
    latexTable = strcat(latexTable, '\textbf{$\lambda_i$} & \textbf{$\Re(\lambda_i)$} & \textbf{$\Im(\lambda_i)$} & \textbf{M. P. St.$_1$} & \textbf{M. P. St.$_2$} & \textbf{M. P. St.$_3$} & \textbf{Nat. Freq. (Hz)} & \textbf{Damp. Freq. (Hz)} & \textbf{Damp. ratio (\%)} \\');
    latexTable = strcat(latexTable, '\hline');

    for i = 1:numel(selected_modes)
        partmode = sort(partMatrix(:,selected_modes(i)),'descend');
        latexTable = strcat(latexTable, sprintf('$%s$ & %0.2f & %0.2f & $%s$ (%0.2f) & $%s$  (%0.2f) & $%s$  (%0.2f) & %0.2f & %0.2f & %0.2f \\', ...
            xlabels{selected_modes(i)}, real(eigenVal(selected_modes(i))), imag(eigenVal(selected_modes(i))), string_values{selected_modes(i),1}, partmode(1)*100, string_values{selected_modes(i),2}, partmode(2)*100, string_values{selected_modes(i),3}, partmode(3)*100,  ...
            undamped_freqs(selected_modes(i)), damped_freqs(selected_modes(i)), damp_ratio(selected_modes(i))*100));
        latexTable = strcat(latexTable, '\\hline');
    end

    latexTable = strcat(latexTable, '\end{tabular}');

    % Create figure and annotation
    figure;
    ha = annotation('textbox',[0.0288,0.5,1,0.2], 'Interpreter', 'latex');
    set(ha, 'String', latexTable,'EdgeColor', 'none', 'FontSize', 14);
end

end