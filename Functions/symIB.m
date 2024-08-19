function [A_IB0,B_IB0,C_IB0,D_IB0,states_IB,inputs_IB,outputs_IB,nStates_IB,nIP_IB,nOP_IB] = symIB()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upstream Network (Infinite Bus)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Small-signal model for upstream network 
% DC side considered 

% Parameters
syms wb Rup Lup

% States (x)
syms theta_up igd_g igq_g

% Inputs (u)
syms wup Vup vgd_g vgq_g

%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

stateVec         = [igd_g igq_g theta_up];
us_vec           = [wup Vup];
ug_vec           = [vgd_g vgq_g];
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


%%%%%%%%%%%%%%%% State equations %%%%%%%%%%%%%%%
% First converting input voltages from global reference frame to local one
% - Assumimg that upstream network will be the slack node if ever used:
vgd         = vgd_g; 
vgq         = vgq_g;


digd_g     = (wb/Lup)*(Vup - vgd - Rup*igd_g + wup*Lup*igq_g);
digq_g     = (wb/Lup)*(- vgq - Rup*igq_g - wup*Lup*igd_g);

dtheta_up   = wb*wup;

diffeqVec   = [digd_g digq_g dtheta_up];

%%%%%%%%%%%%%%%% Output vector %%%%%%%%%%%%%%%

p_up              =  Vup*igd_g;
q_up              = -Vup*igq_g;

% Output currents from local referemce frame to global one
% Assumimg that upstream network will be the slack node if ever used:

outputeqVec_g = [igd_g igq_g theta_up wup];
outputeqVec_s = [p_up q_up];

outputeqVec      = [outputeqVec_s outputeqVec_g];

Fx          = simplify(jacobian(diffeqVec, stateVec));
Fu          = simplify(jacobian(diffeqVec, inputVec));

Hx          = simplify(jacobian(outputeqVec, stateVec));
Hu          = simplify(jacobian(outputeqVec, inputVec));

Ai          = Fx ;
Bi          = Fu ;
Ci          = Hx ;
Di          = Hu ;

A_IB0       = subs(Ai, [stateVec inputVec ], [stateVecEq inputVecEq ]);
B_IB0       = subs(Bi, [stateVec inputVec ], [stateVecEq inputVecEq ]);
C_IB0       = subs(Ci, [stateVec inputVec ], [stateVecEq inputVecEq ]);
D_IB0       = subs(Di, [stateVec inputVec ], [stateVecEq inputVecEq ]);

states_IB   = ["i_{d" "i_{q" "\theta_{" ]; % the states names
inputs_IB   = ["\omega" "V_{up}"]; % the local inputs names
outputs_IB   = ["p_e" "q_e"]; % the local outputs names
nStates_IB  = size(A_IB0,1); 
nIP_IB      = [size(B_IB0,2) length(us_vec) length(ug_vec)]; 
nOP_IB      = [size(C_IB0,1) length(outputeqVec_s) length(outputeqVec_g)];  

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssIB_matrices.mat'),'A_IB0','B_IB0','C_IB0','D_IB0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_IB0./A_IB0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_IB0.txt'),'wt');
for i = 1:length(A_IB0)
    for j = 1:length(A_IB0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_IB0(i,j)));            
        end
    end
end
fclose(fid);

end