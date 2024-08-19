function [nodes_Un,nodes_Type,GFM_nType,SG_nType,TR_Un1,TR_Un2] = nodes_NP (Y_network1,Y_DER,Y_TR,Lines_BP,IB_BP,SM_BP,VSC_BP,n_GFM,n_SG)

nodes_Un = zeros(size(Y_network1,1)+size(Y_DER,1),1);
nodes_Type = strings(size(Y_network1,1)+size(Y_DER,1),1);
GFM_nType = strings(n_GFM,1);
SG_nType = strings(n_SG,1);
n_nodesR = size(Y_network1,1);
n_DER = size(Y_DER,1);
Lines_Un = Lines_BP(1);
IB_Un    = IB_BP(1);
SM_Un    = SM_BP(1);
VSC_Un   = VSC_BP(1);
n_TR = size(Y_TR,1);
TR_Un1 = ones(n_TR,1);
TR_Un2 = ones(n_TR,1);

nodes_Un (1:n_nodesR) = Lines_Un;
nodes_Type (1:n_nodesR) = "PQ";

GFM = 0;
SG = 0;

for i = 1:n_DER
    if Y_DER(i,2) == 1      %it is an IB
        nodes_Un (i+n_nodesR) = IB_Un;
        nodes_Type (i+n_nodesR) = "swing";
    elseif Y_DER(i,2) == 2  %it is an GFM
        GFM = GFM + 1;
        nodes_Un (i+n_nodesR) = VSC_Un;
        if Y_DER(i,4) == 1 
            nodes_Type (i+n_nodesR) = "swing";
            GFM_nType(GFM) = "swing"; 
        elseif Y_DER(i,4) == 2
            nodes_Type (i+n_nodesR) = "PV";
            GFM_nType(GFM) = "PV";
        end
    elseif Y_DER(i,2) == 3  %it is an GFL 
        nodes_Un (i+n_nodesR) = VSC_Un;
        nodes_Type (i+n_nodesR) = "PQ";
    elseif Y_DER(i,2) == 4  %it is an SM 
        SG = SG + 1;
        nodes_Un (i+n_nodesR) = SM_Un;
        if Y_DER(i,4) == 1 
            nodes_Type (i+n_nodesR) = "swing";
            SG_nType(SG) = "swing"; 
        elseif Y_DER(i,4) == 2
            nodes_Type (i+n_nodesR) = "PV";
            SG_nType(SG) = "PV";
        end        
    end
end

for i = 1:n_TR
    TR_Un1(i) = max(nodes_Un(Y_TR(i,2)),nodes_Un(Y_TR(i,3)));
    TR_Un2(i) = min(nodes_Un(Y_TR(i,2)),nodes_Un(Y_TR(i,3)));
end

end