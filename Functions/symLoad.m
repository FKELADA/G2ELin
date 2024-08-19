function [A_Load0,B_Load0,C_Load0,D_Load0,states_Load,inputs_Load,outputs_Load,nStates_Load,nIP_Load,nOP_Load] = symLoad()

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(strcat(cd,'/Symbolic'));

% Small-signal model for a RL-series Load branch 

% Parameters
syms wb Rc Lc

% States (x)
syms icd_g icq_g

% %Algebraic (z)
% syms md mq w_pll

% Inputs (u)
syms wg vgd_g vgq_g                  %until vdc_ref --> u_k^s (closed loop inputs), after: u^g


%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

stateVec         = [icd_g icq_g];
us_vec           = [];
ug_vec           = [wg vgd_g vgq_g];
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

dicd_g        = (wb/Lc)*(vgd_g - Rc*icd_g + wg*Lc*icq_g);
dicq_g        = (wb/Lc)*(vgq_g - Rc*icq_g - wg*Lc*icd_g);

diffeqVec   = [dicd_g dicq_g];

%%%%%%%%%%%%%%%%% Output Vector: %%%%%%%%%%%%%%%%%%%%%%
outputeqVec_g = [icd_g icq_g];

p_c              = (vgd_g*icd_g + vgq_g*icq_g);
q_c              = (-vgd_g*icq_g + vgq_g*icd_g);

outputeqVec_s = [];

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

A_Load0      = subs(Ai, [stateVec inputVec ], [stateVecEq inputVecEq]);
B_Load0      = subs(Bi, [stateVec inputVec ], [stateVecEq inputVecEq]);
C_Load0      = subs(Ci, [stateVec inputVec ], [stateVecEq inputVecEq]);
D_Load0      = subs(Di, [stateVec inputVec ], [stateVecEq inputVecEq]);

states_Load  = ["i_{{c_d}" "i_{{c_q}"];
inputs_Load  = ""; % The local inputs names
% outputs_Load = ["i_{c_d}" "i_{c_q}" "p_c" "q_c"]; % The local outputs names
outputs_Load = ""; % The local outputs names
nStates_Load = size(A_Load0,1); 
nIP_Load     = [size(B_Load0,2) length(us_vec) length(ug_vec)]; 
nOP_Load     = [size(C_Load0,1) length(outputeqVec_s) length(outputeqVec_g)];

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssLoad_matrices.mat'),'A_Load0','B_Load0','C_Load0','D_Load0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_Load0./A_Load0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_Load0.txt'),'wt');
for i = 1:length(A_Load0)
    for j = 1:length(A_Load0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_Load0(i,j)));            
        end
    end
end
fclose(fid);

end