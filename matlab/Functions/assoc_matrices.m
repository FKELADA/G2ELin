function [F,G,K,L] = assoc_matrices (N_s,N_G,N_CF,N_Cf,Y_network,Y_line1,n_total_inputs,n_total_inputs_s,n_total_inputs_g,n_total_outputs,n_total_outputs_s,n_total_outputs_g,n_nodes,n_DG,n_GFM,n_GFL,n_SG,nIP_s,nIP_GFL,nIP_GFM,nIP_SG,nIP_Node,nIP_Line,nIP_Load,sflag_GFM,sflag_SG,nOP_s,nOP_GFM,nOP_GFL,nOP_SG,nOP_Line,nOP_Load,nOP_Node,n_lines,n_loads)
%% The F_Matrix:
% The F_matrix maps the algebraic relations between all the concatenated
% inputs of all state spaces and all the concatented ((local)) inputs of all state
% spaces as well.
F = zeros(n_total_inputs,n_total_inputs_s);

% For the slack DG:
F(1:nIP_s(2),1:nIP_s(2)) = eye(nIP_s(2));

% For all other DGs:
row_index = nIP_s(2)+2;
col_index = nIP_s(2);
for i = 2:n_DG
    if Y_network((n_nodes-n_DG)+i,2) == 2  %it is a GFM  
        F(row_index+1:row_index+nIP_GFM(2),col_index+1:col_index+nIP_GFM(2)) = eye(nIP_GFM(2));
        row_index = row_index + nIP_GFM(2) + 3;
        col_index = col_index + nIP_GFM(2);
    elseif Y_network((n_nodes-n_DG)+i,2) == 3 %it is a GFL
        F(row_index+1:row_index+nIP_GFL(2),col_index+1:col_index+nIP_GFL(2)) = eye(nIP_GFL(2));
        row_index = row_index + nIP_GFL(2) + 3;
        col_index = col_index + nIP_GFL(2);
    elseif Y_network((n_nodes-n_DG)+i,2) == 4 %it is a SG
        F(row_index+1:row_index+nIP_SG(2),col_index+1:col_index+nIP_SG(2)) = eye(nIP_SG(2));
        row_index = row_index + nIP_SG(2) + 3;
        col_index = col_index + nIP_SG(2);
    end
end

% figure
% heatmap(F)
% the end of the F_matrix
%% The G_matrix:
% The G_matrix maps the algebraic relations between all the concatenated
% inputs of all state spaces and all the concatented outputs of all state
% spaces as well.
G = zeros(n_total_inputs,n_total_outputs);
% General preparations:
% Notice the order is: slack -> DGs -> nodes -> lines -> loads
% start & ending rows of slack/DG/Nodes/lines/loads:
               %[start   %end]
rows_slack    = [1   nIP_s(1)];
col_slack     = [1   nOP_s(1)];
rows_DG       = [rows_slack(2)+1   rows_slack(2)+nIP_SG(1)*(n_SG-sflag_SG)+nIP_GFM(1)*(n_GFM-sflag_GFM)+nIP_GFL(1)*n_GFL];
col_DG        = [col_slack(2)+1    col_slack(2)+nOP_SG(1)*(n_SG-sflag_SG)+nOP_GFM(1)*(n_GFM-sflag_GFM)+nOP_GFL(1)*n_GFL];
rows_nodes    = [rows_DG(2)+1    rows_DG(2)+(n_nodes-n_DG)*nIP_Node(1)];
col_nodes     = [col_DG(2)+1     col_DG(2)+(n_nodes-n_DG)*nOP_Node(1)];
rows_Lines    = [rows_nodes(2)+1    rows_nodes(2)+n_lines*nIP_Line(1)];
col_Lines     = [col_nodes(2)+1     col_nodes(2)+n_lines*nOP_Line(1)];
rows_Loads    = [rows_Lines(2)+1    rows_Lines(2)+n_loads*nIP_Load(1)];
col_Loads     = [col_Lines(2)+1     col_Lines(2)+n_loads*nOP_Load(1)];

% thus slack rows start at rows_slack(1) and ends at rows_slack(2) and
% columns start at col_slack(1) and ends at col_slack(2) and so on ..

%%%%% Mapping relations of the slack DG:
% And Since the slack node is always the first one:
% The only subsitution is for its network node voltage thus:
% Part I:
G(nIP_s(2)+1:nIP_s(2)+2,col_DG(2)+2*(N_s-1)+1:col_DG(2)+2*(N_s-1)+2) = eye(2);

%%%%% Mapping relations of all the other DGs (by order of appearance in the network):
% N.B.: an Infinite bus (upstream network) can never be a normal DG, it is
% always assumed to be the slack DG AND always at node 1 of a given network
% topology if ever used.

% Part II:
% first: Mapping the thetas relations:
% remember all other DGs other than slack one have an input theta that
% comes from the slack one (an output of the slack state space), thus:
col_index = nOP_s(2)+3;
row_index = nIP_s(2);
for i  = 2:n_DG
    if Y_network((n_nodes-n_DG)+i,2) == 2  %it is a GFM  
        row_index = row_index + nIP_GFM(2) + 3;
    elseif Y_network((n_nodes-n_DG)+i,2) == 3 %it is a GFL
        row_index = row_index + nIP_GFL(2) + 3;
    elseif Y_network((n_nodes-n_DG)+i,2) == 4 %it is a SG
        row_index = row_index + nIP_SG(2) + 3; 
    end
    G(row_index,col_index) = 1;
end

% Part III:
% Second: Mapping the DG node voltages relations:
SG = sflag_SG; GFM = sflag_GFM; GFL = 0;
row_index = rows_slack(2)-2;
for i  = 2:n_DG
    if Y_network((n_nodes-n_DG)+i,2) == 2  %it is a GFM  
        GFM = GFM + 1;
        row_index = row_index + nIP_GFM(2) + 3;
        col_index = col_DG(2) + 2*(N_CF(GFM)-1) + 1;
    elseif Y_network((n_nodes-n_DG)+i,2) == 3 %it is a GFL
        GFL = GFL + 1;
        row_index = row_index + nIP_GFL(2) + 3;
        col_index = col_DG(2) + 2*(N_Cf(GFL)-1) + 1;
    elseif Y_network((n_nodes-n_DG)+i,2) == 4 %it is a SG
        SG = SG + 1;
        row_index = row_index + nIP_SG(2) + 3;
        col_index = col_DG(2) + 2*(N_G(SG)-1) + 1;
    end
    G(row_index+1:row_index+2,col_index:col_index+1) = eye(2);
end

%%%%% Maping relations of the Nodes state spaces:
% remember nodes' currents are algebraically related to (1) currents from DGs, (2) currents from
% lines and (3) currents from loads states spaces. These relations are
% described by Part V of this code. However Part IV, descibes the wg input
% from the slack DG to all nodes' state spaces.

% Part IV:
% Mapping the omega relations:
col_index = col_slack(2);
row_index = rows_nodes(1); 
for i = 1:n_nodes-n_DG
    G(row_index,col_index) = 1;
    row_index = row_index + 3;
end

% Part V:
% Mapping the nodes currents relations:
% Convention: sum Iin = sum I out -> thus always Ig from DGs are always in
% to the node and Ic from loads are laways out of the node -> Ig always +1
% and Ic always -1 if any
% However it is a little bit trickier for lines: if a line is defined in
% the Y_line matrix as "from" a certain node: it goes out of that node thus
% considered as -1 and if it is defined as "to" from the specific node it
% will be +1:

% first for the DG currents (if any):
% It should be noted that all DGs are declared as the first (n_DG) nodes in
% the Y_network matrix, thus:
% first for the slack then for all other DGs:
G(rows_DG(2)+3*(N_s-1)+2:rows_DG(2)+3*(N_s-1)+3,nOP_s(2)+1:nOP_s(2)+2) = eye(2);

% for all other DGs:
col_index = col_slack(2)-2;
i_GFM = sflag_GFM; i_SG = sflag_SG; i_GFL = 0;
for i = 2:n_DG 
    if Y_network((n_nodes-n_DG)+i,2) == 2  %it is an GFM
        i_GFM = i_GFM + 1;
        row_index = rows_nodes(1) + 3*(N_CF(i_GFM)-1);
        col_index = col_index + nOP_GFM(2) + 2;
    elseif Y_network((n_nodes-n_DG)+i,2) == 3  %it is an GFL
        i_GFL = i_GFL + 1;
        row_index = rows_nodes(1) + 3*(N_Cf(i_GFL)-1);       
        col_index = col_index + nOP_GFL(2) + 2;
    elseif Y_network((n_nodes-n_DG)+i,2) == 4  %it is an SG
        i_SG = i_SG + 1;
        row_index = rows_nodes(1) + 3*(N_G(i_SG)-1);        
        col_index = col_index + nOP_SG(2) + 2;
    end
    G(row_index+1:row_index+2,col_index+1:col_index+2) = eye(2);
end

% second for the lines currents (+1 if "to" the node and -1 if "from" the node):
row_index = rows_nodes(1); 
for i = 1:n_nodes-n_DG
    [indices_f flag_f] = find((Y_line1(:,2)) == i); % find the number of lines coming "from" this node and their indices (1st line, 2nd line or ..etc)
    nlines_f = length(flag_f);
    [indices_t flag_t] = find((Y_line1(:,3)) == i); % find the number of lines going "to" this node and their indices (1st line, 2nd line or ..etc)
    nlines_t = length(flag_t);    
    for j = 1:nlines_f
        col_index = col_nodes(2) - 2 + 2*indices_f(j); %indices_f(j): is the index of the jth line
        G(row_index+1:row_index+2,col_index+1:col_index+2) = -eye(2); 
    end
    for k = 1:nlines_t
        col_index = col_nodes(2) - 2 + 2*indices_t(k); %indices_t(k): is the index of the kth line
        G(row_index+1:row_index+2,col_index+1:col_index+2) = eye(2); 
    end   
    row_index = row_index + 3;
end

% third for the load currents (if any):
%N.B.: In declaring loads the following convention must be used: if Load 1
%is at node i thus, load 2 must be at node i+n where n>0 
row_index = rows_nodes(1); 
Load = 0;
for i = 1:n_nodes-n_DG
    if Y_network(i,3) == 1 % find if their is a load connected to thins node and its index (1st load, 2nd load or ..etc)
        Load = Load + 1;
        col_index = col_Lines(2) - 2 + 2*Load; %Load is the index of the load
        G(row_index+1:row_index+2,col_index+1:col_index+2) = -eye(2); 
    end 
    row_index = row_index + 3;
end

% Part VI:
% Mapping the lines' voltages to the respective nodes voltages:
% Remember lines have respective nodes' voltages as inputs to their state
% spaces. This is where these relations are defined.

% first: mapping the input omega_g from slack nodes to the lines' state
% spaces:
col_index = col_slack(2);
row_index = rows_Lines(1); 
for i = 1:n_lines
    G(row_index,col_index) = 1;
    row_index = row_index + 5;
end

% second: mapping the lines' voltages to nodes' voltages:
row_index = rows_Lines(1); 
for i = 1:n_lines
    col_index = col_DG(2) - 2 + 2*(Y_line1(i,2)); %(Y_line(i,2)-n_DG) is the node where the line goes "from" and (Y_line(i,3)-n_DG) is the node where the line goes "to"
    G(row_index+1:row_index+2,col_index+1:col_index+2) = eye(2);
    row_index = row_index + 5;
end


row_index = rows_Lines(1)+2; 
for i = 1:n_lines
    col_index = col_DG(2) - 2 + 2*(Y_line1(i,3)); %(Y_line(i,2)-n_DG) is the node where the line goes "from" and (Y_line(i,3)-n_DG) is the node where the line goes "to"
    G(row_index+1:row_index+2,col_index+1:col_index+2) = eye(2);
    row_index = row_index + 5;
end

% Part VII:
% Mapping the loads' voltages to the respective nodes voltages:
% Remember loads have respective nodes' voltages as inputs to their state
% spaces. This is where these relations are defined.

% first: mapping the input omega_g from slack nodes to the loads' state
% spaces:
col_index = col_slack(2);
row_index = rows_Loads(1); 
for i = 1:n_loads
    G(row_index,col_index) = 1;
    row_index = row_index + 3;
end

% second: mapping the lines' voltages to nodes' voltages:
row_index = rows_Loads(1); 
for ii = 1:n_nodes-n_DG
    if Y_network(ii,3) == 1 
        col_index = col_DG(2) - 2 + 2*ii; %ii is the index of the node
        G(row_index+1:row_index+2,col_index+1:col_index+2) = eye(2); 
        row_index = row_index + 3;
    end 
end

% figure
% heatmap(G)
% the end of the G matrix

%% The K_matrix:
K = zeros(n_total_outputs_s,n_total_inputs_s);

% figure
% heatmap(K)
% the end of the K_matrix
%% The L_matrix:
% The L_matrix maps the algebraic relations between all the concatenated
% ((Local)) outputs of all state spaces and all the concatented outputs of all state
% spaces as well.
L = zeros(n_total_outputs_s,n_total_outputs);

% For the slack DG:
L(1:nOP_s(2),1:nOP_s(2)) = eye(nOP_s(2));

% For all other DGs:
row_index = nOP_s(2);
col_index = nOP_s(2)+4;
for i = 2:n_DG
    if Y_network((n_nodes-n_DG)+i,2) == 2  %it is a GFM  
        L(row_index+1:row_index+nOP_GFM(2),col_index+1:col_index+nOP_GFM(2)) = eye(nOP_GFM(2));
        row_index = row_index + nOP_GFM(2);
        col_index = col_index + nOP_GFM(2) + 2;
    elseif Y_network((n_nodes-n_DG)+i,2) == 3 %it is a GFL
        L(row_index+1:row_index+nOP_GFL(2),col_index+1:col_index+nOP_GFL(2)) = eye(nOP_GFL(2));
        row_index = row_index + nOP_GFL(2);
        col_index = col_index + nOP_GFL(2) + 2;
    elseif Y_network((n_nodes-n_DG)+i,2) == 4 %it is a SG
        L(row_index+1:row_index+nOP_SG(2),col_index+1:col_index+nOP_SG(2)) = eye(nOP_SG(2));
        row_index = row_index + nOP_SG(2);
        col_index = col_index + nOP_SG(2) + 2;
    end
end

% figure
% heatmap(L)
% the end of the L_matrix
%%
% ylabels = inputNames;
% xlabels = outputNames;
% figure
% heatmap(xlabels,ylabels,G)
% G11      = G(rows_slack(1):rows_DG(2),col_slack(1):col_DG(2));
% figure
% heatmap(G11)
% G12      = G(rows_slack(1):rows_DG(2),col_DG(2)+1:col_nodes(2));
% figure
% heatmap(G12)
% G13      = G(rows_slack(1):rows_DG(2),col_nodes(2)+1:col_Lines(2));
% figure
% heatmap(G13)
% G13      = G(rows_slack(1):rows_DG(2),col_nodes(2)+1:col_Lines(2));
% figure
% heatmap(G13)
% heatmap(G13)
% G14      = G(rows_slack(1):rows_DG(2),col_Lines(2)+1:col_Loads(2));
% figure
% heatmap(G14)
% G13      = G(1:n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line,n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line + 1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line + n_loads*nOP_Load);
% G21      = G(n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + 1 : n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line,1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line);
% G22      = G(n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + 1 : n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line,n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + 1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line);
% G23      = G(n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + 1 : n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line,n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line + 1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line + n_loads*nOP_Load);
% G31      = G(n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line + 1 : n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line + n_loads*nIP_Load,1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line);
% G32      = G(n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line + 1 : n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line + n_loads*nIP_Load,n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + 1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line);
% G33      = G(n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line + 1 : n_total_inputs - n_loads*nIP_Load - n_lines*nIP_Line + n_lines*nIP_Line + n_loads*nIP_Load,n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line + 1:n_total_outputs - n_loads*nOP_Load - n_lines*nOP_Line + n_lines*nOP_Line + n_loads*nOP_Load);

end