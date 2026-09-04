function [F,G,K,L] = assoc_matrices_s (N_s,N_G,N_CF,N_Cf,Y_network,Y_line1,Y_bus,n_total_inputs,n_total_inputs_s,n_total_inputs_g,n_total_outputs,n_total_outputs_s,n_total_outputs_g,n_nodes,n_DG,n_GFM,n_GFL,n_SG,nIP_s,nIP_GFL,nIP_GFM,nIP_SG,nIP_Node,nIP_Line,nIP_Load,sflag_GFM,sflag_SG,nOP_s,nOP_GFM,nOP_GFL,nOP_SG,nOP_Line,nOP_Load,nOP_Node,n_lines,n_loads)
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

% thus slack rows start at rows_slack(1) and ends at rows_slack(2) and
% columns start at col_slack(1) and ends at col_slack(2) and so on ..

%%%%% Mapping relations of the slack DG:
% And Since the slack node is always the first one:
% The only subsitution is for its network node voltage thus:
% Part I:
G(nIP_s(2)+1:nIP_s(2)+2,nOP_s(2)+2*(N_s-1)+1:nOP_s(2)+2*(N_s-1)+2) = inv(abs(Y_bus(N_s,N_s)));

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
col_index = col_slack(2);
for i  = 2:n_DG
    if Y_network((n_nodes-n_DG)+i,2) == 2  %it is a GFM  
        GFM = GFM + 1;
        row_index = row_index + nIP_GFM(2) + 3;
        col_index = col_index + nOP_GFM(2) + 1;
        G(row_index+1:row_index+2,col_index:col_index+1) = inv(abs(Y_bus(N_CF(GFM),N_CF(GFM))));
    elseif Y_network((n_nodes-n_DG)+i,2) == 3 %it is a GFL
        GFL = GFL + 1;
        row_index = row_index + nIP_GFL(2) + 3;
        col_index = col_index + nOP_GFL(2) + 1;
        G(row_index+1:row_index+2,col_index:col_index+1) = inv(abs(Y_bus(N_Cf(GFL),N_Cf(GFL))));
    elseif Y_network((n_nodes-n_DG)+i,2) == 4 %it is a SG
        SG = SG + 1;
        row_index = row_index + nIP_SG(2) + 3;
        col_index = col_index + nOP_SG(2) + 1;
        G(row_index+1:row_index+2,col_index:col_index+1) = inv(abs(Y_bus(N_G(SG),N_G(SG))));
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