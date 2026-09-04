function [N_s, N_G, N_CF, N_Cf, N_L] = nodes_set(Y_network,Y_DER)
% N_s is the node at which the slack DG is connected
% N_G are the nodes at which a SM is connected
% N_CF are the nodes at which a GFM is connected
% N_Cf are the nodes at which a GFL is connected
% N_L are the nodes at which a load is connected

N_s = [];
N_G = [];
N_CF = [];
N_Cf = [];

n_sg = 0;
n_gfm = 0;
n_gfl = 0;

for i = 1:size(Y_DER,1)
    if Y_DER(i,4) == 1 %it's type is slack (it is only 1)
        N_s = Y_DER(i,1);
    end
end

for i = 1:size(Y_DER,1)
    if Y_DER(i,2) == 2 %it is a GFM
        n_gfm = n_gfm + 1;
        N_CF(n_gfm) = Y_DER(i,1);
    elseif Y_DER(i,2) == 3 %it is a GFL
        n_gfl = n_gfl + 1;
        N_Cf(n_gfl) = Y_DER(i,1);
    elseif Y_DER(i,2) == 4 %it is a SM
        n_sg = n_sg + 1;
        N_G(n_sg) = Y_DER(i,1);
    end
end

N_L = [];
n_load = 0;
for i = 1:size(Y_network,1)
    if Y_network(i,3) == 1 %there is a load connected to this node
        n_load = n_load + 1;
        N_L(n_load) = Y_network(i,1);
    end
end

end