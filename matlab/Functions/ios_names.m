function [inputNames,outputNames,stateNames] = ios_names(n_nodes,n_DG,inputs_s,inputs_GFL,inputs_GFM,inputs_SG,outputs_s,outputs_GFL,outputs_GFM,outputs_SG,outputs_Node,outputs_Line,outputs_Load,n_total_inputs_s,n_total_outputs_s,n_GFM,n_GFL,n_SG,nIP_GFL,nIP_GFM,nIP_SG,nIP_Node,sflag_GFM,sflag_SG,nOP_GFM,nOP_GFL,nOP_SG,nOP_Line,nOP_Load,nOP_Node,states_s,states_GFL,states_GFM,states_SG,n_lines,n_loads,states_Line,states_Load,states_Node,nStates_s,nStates_SG,nStates_GFM,nStates_GFL,nStates_Line,nStates_Load,nStates_Node)

    inputNames_s  = strcat(inputs_s,'_s');
    outputNames_s = strcat(outputs_s,'_s');
    for i = 1:n_GFL
        inputNames_nGFL1(i,:) = strcat(strtrim(cellstr(inputs_GFL)),'_{GFL_',strtrim(cellstr(num2str(i))),'}');
        outputNames_nGFL1(i,:) = strcat(strtrim(cellstr(outputs_GFL)),'_{GFL_',strtrim(cellstr(num2str(i))),'}');
    end 
    for i = 1:(n_GFM-sflag_GFM)
        inputNames_nGFM1(i,:) = strcat(strtrim(cellstr(inputs_GFM)),'_{GFM_',strtrim(cellstr(num2str(i))),'}');
        outputNames_nGFM1(i,:) = strcat(strtrim(cellstr(outputs_GFM)),'_{GFM_',strtrim(cellstr(num2str(i))),'}');
    end
    for i = 1:(n_SG-sflag_SG)
        inputNames_nSG1(i,:) = strcat(strtrim(cellstr(inputs_SG)),'_{SG_',strtrim(cellstr(num2str(i))),'}');
        outputNames_nSG1(i,:) = strcat(strtrim(cellstr(outputs_SG)),'_{SG_',strtrim(cellstr(num2str(i))),'}');
    end

    for i = 1:(n_nodes-n_DG)
        outputNames_nNode1(i,:) = strcat(strtrim(cellstr(outputs_Node)),'_{Nd_',strtrim(cellstr(num2str(i))),'}');
    end
    for i = 1:n_loads
        outputNames_nLoad1(i,:) = strcat(strtrim(cellstr(outputs_Load)),'_{Ld_',strtrim(cellstr(num2str(i))),'}');
    end
    for i = 1:n_lines
        outputNames_nLine1(i,:) = strcat(strtrim(cellstr(outputs_Line)),'_{Ln_',strtrim(cellstr(num2str(i))),'}');
    end
    inputNames = strings(1,n_total_inputs_s);
    outputNames = strings(1,n_total_outputs_s);
    
    ii=1;
    if n_GFL > 0
    for j =1:nIP_GFL(2):nIP_GFL(2)*n_GFL
        inputNames_nGFL2(1,j:j+(nIP_GFL(2)-1)) = inputNames_nGFL1(ii,:) ;
        ii = ii + 1;
    end
    else inputNames_nGFL2 = strings;
    end
    ii=1;
    if n_GFL > 0
    for j =1:nOP_GFL(2):nOP_GFL(2)*n_GFL
        outputNames_nGFL2(1,j:j+(nOP_GFL(2)-1)) = outputNames_nGFL1(ii,:) ;
        ii = ii + 1;
    end
    else outputNames_nGFL2 = strings;
    end
    ii=1;
    if (n_GFM-sflag_GFM) > 0
    for j =1:nIP_GFM(2):nIP_GFM(2)*(n_GFM-sflag_GFM)
        inputNames_nGFM2(1,j:j+(nIP_GFM(2)-1)) = inputNames_nGFM1(ii,:) ;
        ii = ii + 1;
    end
    else inputNames_nGFM2 = strings;
    end
    ii=1;
    if (n_GFM-sflag_GFM) > 0
    for j =1:nOP_GFM(2):nOP_GFM(2)*(n_GFM-sflag_GFM)
        outputNames_nGFM2(1,j:j+(nOP_GFM(2)-1)) = outputNames_nGFM1(ii,:) ;
        ii = ii + 1;
    end
    else outputNames_nGFM2 = strings;
    end
    ii=1;
    if (n_SG-sflag_SG) > 0
    for j =1:nIP_SG(2):nIP_SG(2)*(n_SG-sflag_SG)
        inputNames_nSG2(1,j:j+(nIP_SG(2)-1)) = inputNames_nSG1(ii,:) ;
        ii = ii + 1;
    end
    else inputNames_nSG2 = strings;
    end
    ii=1;
    if (n_SG-sflag_SG) > 0
    for j =1:nOP_SG(2):nOP_SG(2)*(n_SG-sflag_SG)
        outputNames_nSG2(1,j:j+(nOP_SG(2)-1)) = outputNames_nSG1(ii,:) ;
        ii = ii + 1;
    end
    else outputNames_nSG2 = strings;
    end
    
    ii=1;
    for j =1:nOP_Node(2):nOP_Node(2)*n_nodes
        outputNames_nNode2(1,j:j+(nOP_Node(2)-1)) = outputNames_nNode1(ii,:) ;
        ii = ii + 1;
    end
    ii=1;
    for j =1:nOP_Line(2):nOP_Line(2)*n_lines
        outputNames_nLine2(1,j:j+(nOP_Line(2)-1)) = outputNames_nLine1(ii,:) ;
        ii = ii + 1;
    end
    ii=1;
    for j =1:nOP_Load(2):nOP_Load(2)*n_loads
        outputNames_nLoad2(1,j:j+(nOP_Load(2)-1)) = outputNames_nLoad1(ii,:) ;
        ii = ii + 1;
    end
    inputNames = [inputNames_s inputNames_nGFM2 inputNames_nSG2 inputNames_nGFL2];
    outputNames = [outputNames_s outputNames_nGFM2 outputNames_nSG2 outputNames_nGFL2];
    % Find empty positions using the "isempty" function
    emptyPositionsIP = cellfun(@isempty, inputNames);
    emptyPositionsOP = cellfun(@isempty, outputNames);
    % Remove empty positions using logical indexing
    inputNames = inputNames(~emptyPositionsIP);
    outputNames = outputNames(~emptyPositionsOP);
    
    stateNames_s = strcat(states_s,'_s}');
    for i = 1:n_GFL
        stateNames_nGFL1(i,:) = strcat(strtrim(cellstr(states_GFL)),'_{GFL_',strtrim(cellstr(num2str(i))),'}}');
    end 
    for i = 1:(n_GFM-sflag_GFM)
        stateNames_nGFM1(i,:) = strcat(strtrim(cellstr(states_GFM)),'_{GFM_',strtrim(cellstr(num2str(i))),'}}');
    end
    for i = 1:(n_SG-sflag_SG)
        stateNames_nSG1(i,:) = strcat(strtrim(cellstr(states_SG)),'_{SG_',strtrim(cellstr(num2str(i))),'}}');
    end
    for i = 1:(n_nodes-n_DG)
        stateNames_nNode1(i,:) = strcat(strtrim(cellstr(states_Node)),'_{{Nd}_',strtrim(cellstr(num2str(i))),'}}');
    end     
    for i = 1:n_lines
        stateNames_nLine1(i,:) = strcat(strtrim(cellstr(states_Line)),'_{{Ln}_',strtrim(cellstr(num2str(i))),'}}');
    end 
    for i = 1:n_loads
        stateNames_nLoad1(i,:) = strcat(strtrim(cellstr(states_Load)),'_{{Ld}_',strtrim(cellstr(num2str(i))),'}}');
    end 
    stateNames = strings(1,nStates_s+nStates_GFM*(n_GFM-sflag_GFM)+nStates_SG*(n_SG-sflag_SG)+nStates_GFL*n_GFL+nStates_Line*n_lines+nStates_Load*n_loads);
    ii=1;
    if n_GFL > 0
    for j =1:nStates_GFL:nStates_GFL*n_GFL
        stateNames_nGFL2(1,j:j+(nStates_GFL-1)) = stateNames_nGFL1(ii,:) ;
        ii = ii + 1;
    end
    else stateNames_nGFL2 = strings;
    end
    ii=1;
    if (n_GFM-sflag_GFM) > 0
    for j =1:nStates_GFM:nStates_GFM*(n_GFM-sflag_GFM)
        stateNames_nGFM2(1,j:j+(nStates_GFM-1)) = stateNames_nGFM1(ii,:) ;
        ii = ii + 1;
    end
    else stateNames_nGFM2 = strings;
    end
    ii=1;
    if (n_SG-sflag_SG) > 0
    for j =1:nStates_SG:nStates_SG*(n_SG-sflag_SG)
        stateNames_nSG2(1,j:j+(nStates_SG-1)) = stateNames_nSG1(ii,:) ;
        ii = ii + 1;
    end
    else stateNames_nSG2 = strings;
    end
    ii=1;
    for j =1:nStates_Line:nStates_Line*n_lines
        stateNames_nLine2(1,j:j+(nStates_Line-1)) = stateNames_nLine1(ii,:) ;
        ii = ii + 1;
    end
    ii=1;
    for j =1:nStates_Load:nStates_Load*n_loads
        stateNames_nLoad2(1,j:j+(nStates_Load-1)) = stateNames_nLoad1(ii,:) ;
        ii = ii + 1;
    end
    ii=1;
    for j =1:nStates_Node:nStates_Node*(n_nodes-n_DG)
        stateNames_nNode2(1,j:j+(nStates_Node-1)) = stateNames_nNode1(ii,:) ;
        ii = ii + 1;
    end
    stateNames = [stateNames_s stateNames_nGFM2 stateNames_nSG2 stateNames_nGFL2 stateNames_nNode2 stateNames_nLine2 stateNames_nLoad2];
    % Find empty positions using the "isempty" function
    emptyPositions = cellfun(@isempty, stateNames);
    % Remove empty positions using logical indexing
    stateNames = stateNames(~emptyPositions);

end