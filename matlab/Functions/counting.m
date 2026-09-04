function [n_DG,n_IB,n_GFM,n_GFL,n_SG,n_loads,n_nodes,n_lines,n_TR,sflag_IB,sflag_GFM,sflag_SG] = counting(Y_network, Y_line, Y_TR)
n_DG = 0;
for i = 1:size(Y_network,1)
    if Y_network(i,2) ~=0
        n_DG = n_DG + 1;
    else n_DG = n_DG;
    end
end
n_IB  = 0;
n_GFM = 0;
n_GFL = 0;
n_SG  = 0;
for i = 1:size(Y_network,1)
    if Y_network(i,2) ==1  % IB
        n_IB = n_IB + 1;
    elseif Y_network(i,2) ==2  % GFM
        n_GFM = n_GFM + 1;
    elseif Y_network(i,2) ==3  % GFL
        n_GFL = n_GFL + 1;
    elseif Y_network(i,2) ==4  % SG
        n_SG = n_SG + 1;
    end
end
n_loads = 0;
for i = 1:size(Y_network,1)
    if Y_network(i,3) ~=0
        n_loads = n_loads + 1;
    else n_loads = n_loads;
    end
end
n_nodes = size(Y_network,1);
n_TR    = size(Y_TR,1); 
n_lines = size(Y_line,1)-n_TR;


sflag_IB  = 0;
sflag_GFM = 0;
sflag_SG  = 0;
first_DG_index = size(Y_network,1) - n_DG + 1;
if Y_network(first_DG_index,2) == 1   % The slack DG is an IB
    sflag_IB = sflag_IB + 1;
elseif Y_network(first_DG_index,2) == 2   % The slack DG is a GFM
    sflag_GFM = sflag_GFM + 1;
elseif Y_network(first_DG_index,2) == 4   % The slack DG is a SG
    sflag_SG = sflag_SG + 1;
end

end