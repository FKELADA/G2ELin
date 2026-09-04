function [Y_network,Y_line,Y_TR,Y_DG] = network_form(Y_DER,Y_network1,Y_line1,Y_TR1,nodes_from_r,nodes_to_r,GFM_P_control)

Y_network = zeros(size(Y_DER,1)+size(Y_network1,1),10);
Y_TR      = zeros(size(Y_DER,1),7);
n_DER     = size(Y_DER,1);

% Forming the Y_network, Y_line & Y_TR matrices:
%% Forming Y_network:
% just concatenating Y_network1 which is the raw network matrix with the
% Y_DER:

Y_DER1 = Y_DER(:,1:10);
Y_DER1(:,1) = size(Y_network1,1)+1:size(Y_network1,1)+size(Y_DER,1);
Y_network = [Y_network1; Y_DER1];


%% Forming Y_TR:
Y_TR(:,1) = 1:n_DER;
Y_TR(:,4) = Y_DER(:,11);
Y_TR(:,5) = Y_DER(:,12);
Y_TR(:,6) = zeros(n_DER,1);
Y_TR(:,7) = ones(n_DER,1);

[nodes_from,index_f] = intersect(Y_DER(:,1),nodes_from_r.');
[nodes_to,index_t]   = intersect(Y_DER(:,1),nodes_to_r.');

if isempty(nodes_from)
    Y_TR(:,2) = zeros(n_DER,1);
end
for i = 1:length(index_f)
    Y_TR(index_f(i),2) = nodes_from(i);
end
if isempty(nodes_to)
    Y_TR(:,3) = zeros(n_DER,1);
end
for i = 1:length(index_t)
    Y_TR(index_t(i),3) = nodes_to(i);
end

bus_after = size(Y_network1,1)+1:size(Y_network1,1)+n_DER;
for i = 1:n_DER
    if Y_TR(i,2) == 0 
        Y_TR(i,2) = bus_after(i); 
    end
end

for i = 1:n_DER
    if Y_TR(i,3) == 0
        Y_TR(i,3) = bus_after(i);   
    end
end


%% Forming Y_line:;

Y_TR(:,1) = size(Y_line1,1)+1:size(Y_line1,1)+size(Y_DER,1);
Y_line = [Y_line1; Y_TR];

%% Forming Y_DG:
Y_DG = strings(size(Y_network, 1), 4);
Y_DG(:, 1:2) = string(Y_network(:, 1:2));
Y_DG(:, 3:4) = repmat("NA", size(Y_network, 1), 2);

for i = 1:size(Y_network, 1)
    if Y_network(i,2) == 1 %IB bus
        Y_DG(i,3) = "IB";
    elseif Y_network(i,2) == 2 %GFM bus
        Y_DG(i,3) = GFM_P_control;
        Y_DG(i,4) = "Ideal";
    elseif Y_network(i,2) == 3 %GFL bus
        Y_DG(i,3) = "PLL";
        Y_DG(i,4) = "Ideal";
    elseif Y_network(i,2) == 4 %SM bus
        Y_DG(i,3) = "SM";
    end
end


% Y_DG is complementary to Y_network:
%            Node DG-type  subType    DCType   
% Y_DG      = [1    4        "SM"       "NA";];

end