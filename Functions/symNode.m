function [A_Node0,B_Node0,C_Node0,D_Node0,states_Node,inputs_Node,outputs_Node,nStates_Node,nIP_Node,nOP_Node] = symNode()

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(strcat(cd,'/Symbolic'));

% Small-signal model for a network node (line capacitance): 

% Parameters
syms wb Cl

% States (x)
syms vgd_g vgq_g

% %Algebraic (z)
% syms md mq w_pll

% Inputs (u)
syms wg ishd_g ishq_g               


%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

stateVec         = [vgd_g vgq_g];
us_vec           = [];
ug_vec           = [wg ishd_g ishq_g];
inputVec         = [us_vec ug_vec];

% Equilibrium vectors
for i = 1:length(inputVec)
    inputVecEq{i} = sym(strcat(char(inputVec(i)),'_0'));
end
inputVecEq = cell2sym(inputVecEq);

for i = 1:length(stateVec)
    stateVecEq{i} = sym(strcat(char(stateVec(i)),'_0'));
end
stateVecEq = cell2sym(stateVecEq);

% for i = 1:length(algVec)
%     algVecEq{i} = sym(strcat(char(algVec(i)),'_0'));
% end
% algVecEq = cell2sym(algVecEq);


%%%%%%%%%%%%%%%%% The circuit in dq-frame: %%%%%%%%%%%%%%%%%%%%%%

dvgd_g        = (wb/Cl)*(ishd_g + wg*Cl*vgq_g);
dvgq_g        = (wb/Cl)*(ishq_g - wg*Cl*vgd_g);

diffeqVec   = [dvgd_g dvgq_g];

%%%%%%%%%%%%%%%%% Output Vector: %%%%%%%%%%%%%%%%%%%%%%
outputeqVec_s = [];
outputeqVec_g = [vgd_g vgq_g];
outputeqVec = [outputeqVec_s outputeqVec_g];

Fx          = simplify(jacobian(diffeqVec, stateVec));
% Fz          = simplify(jacobian(diffeqVec, algVec));
Fu          = simplify(jacobian(diffeqVec, inputVec));

% Gx          = simplify(jacobian(algeqVec, stateVec));
% Gz          = simplify(jacobian(algeqVec, algVec));
% Gu          = simplify(jacobian(algeqVec, inputVec));
% invGz       = inv(Gz);

Hx          = simplify(jacobian(outputeqVec, stateVec));
% Hz          = simplify(jacobian(outputeqVec, algVec));
Hu          = simplify(jacobian(outputeqVec, inputVec));

Ai          = Fx ;
Bi          = Fu ;
Ci          = Hx ;
Di          = Hu ;

A_Node0      = subs(Ai, [stateVec inputVec ], [stateVecEq inputVecEq]);
B_Node0      = subs(Bi, [stateVec inputVec ], [stateVecEq inputVecEq]);
C_Node0      = subs(Ci, [stateVec inputVec ], [stateVecEq inputVecEq]);
D_Node0      = subs(Di, [stateVec inputVec ], [stateVecEq inputVecEq]);

states_Node  = ["v_{{g_d}" "v_{{g_q}"];
inputs_Node  = ""; % The local inputs names
% outputs_Node = ["v_{g_d}" "v_{g_q}"]; % The local outputs names
outputs_Node = ""; % The local outputs names
nStates_Node = size(A_Node0,1); 
nIP_Node     = [size(B_Node0,2) length(us_vec) length(ug_vec)]; 
nOP_Node     = [size(C_Node0,1) length(outputeqVec_s) length(outputeqVec_g)];

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssNode_matrices.mat'),'A_Node0','B_Node0','C_Node0','D_Node0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_Node0./A_Node0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_Node0.txt'),'wt');
for i = 1:length(A_Node0)
    for j = 1:length(A_Node0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_Node0(i,j)));            
        end
    end
end
fclose(fid);

end