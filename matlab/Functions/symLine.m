function [A_Line0,B_Line0,C_Line0,D_Line0,states_Line,inputs_Line,outputs_Line,nStates_Line,nIP_Line,nOP_Line] = symLine()

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(strcat(cd,'/Symbolic'));

% Small-signal model for a RL-series Line 

% Parameters
syms wb Rl Ll

% States (x)
syms ild_g ilq_g

% %Algebraic (z)
% syms md mq w_pll

% Inputs (u)
syms wg vgdj_g vgqj_g vgdk_g vgqk_g                  %until vdc_ref --> u_k^s (closed loop inputs), after: u^g


%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

stateVec         = [ild_g ilq_g];
us_vec           = [];
ug_vec           = [wg vgdj_g vgqj_g vgdk_g vgqk_g];
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

dild_g        = (wb/Ll)*(vgdj_g - vgdk_g - Rl*ild_g + wg*Ll*ilq_g);
dilq_g        = (wb/Ll)*(vgqj_g - vgqk_g - Rl*ilq_g - wg*Ll*ild_g);

diffeqVec   = [dild_g dilq_g];

%%%%%%%%%%%%%%%%% The output vector: %%%%%%%%%%%%%%%%%%%%%%
outputeqVec_s = [];
outputeqVec_g = [ild_g ilq_g];

outputeqVec = [outputeqVec_s outputeqVec_g];

%%%%%%%%%%%%%%%%% Linearization: %%%%%%%%%%%%%%%%%%%%%%

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

A_Line0      = subs(Ai, [stateVec inputVec ], [stateVecEq inputVecEq]);
B_Line0      = subs(Bi, [stateVec inputVec ], [stateVecEq inputVecEq]);
C_Line0      = subs(Ci, [stateVec inputVec ], [stateVecEq inputVecEq]);
D_Line0      = subs(Di, [stateVec inputVec ], [stateVecEq inputVecEq]);

states_Line  = ["i_{{l_d}" "i_{{l_q}"];
inputs_Line  = ""; % The local inputs names
% outputs_Line = ["i_{l_d}" "i_{l_q}"]; % The local outputs names
outputs_Line = ""; % The local outputs names
nStates_Line = size(A_Line0,1); 
nIP_Line     = [size(B_Line0,2) length(us_vec) length(ug_vec)]; 
nOP_Line     = [size(C_Line0,1) length(outputeqVec_s) length(outputeqVec_g)]; 

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssLine_matrices.mat'),'A_Line0','B_Line0','C_Line0','D_Line0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_Line0./A_Line0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_Line0.txt'),'wt');
for i = 1:length(A_Line0)
    for j = 1:length(A_Line0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_Line0(i,j)));            
        end
    end
end
fclose(fid);

end