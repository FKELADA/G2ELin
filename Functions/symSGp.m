function [A_SG0,B_SG0,C_SG0,D_SG0,states_SG,inputs_SG,outputs_SG,nStates_SG,nIP_SG,nOP_SG] = symSGp(PSS,i_DG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SG #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(strcat(cd,'/Symbolic'));

% PSS_Mode:
PSS_type = PSS;

% Small-signal model for SG+TR connected to PCC 

% Parameters
syms wb Rt Lt Rg Ra Ll Lad Laq Lf1d Lffd Rfd L11d L11q L22q R1d R1q R2q H KD mp TG T_LP T_HP K_PSS T1n T1d T2n T2d Tr Ka Ta Ke Te Kfd Tfd

% States (x)
syms igd igq phi_fd phi_1d phi_1q phi_2q Dwr theta Pm e1 e2 efd e3 Dw1 v1 v2 vpss

%Algebraic (z)
syms wr id iq ved veq phi_d phi_q vgd vgq Cm Ce DP Et

% Inputs (u)
syms v_ref p_ref w_ref theta_g vgd_g vgq_g                  %until vdc_ref --> u_k^s (closed loop inputs), after: u^g

%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

algVec           = [wr id iq ved veq phi_d phi_q vgd vgq Cm Ce DP Et];
us_vec           = [v_ref p_ref w_ref];

if i_DG == 1 
    ug_vec           = [vgd_g vgq_g];
    inputVec         = [us_vec ug_vec];
else 
    ug_vec           = [theta_g vgd_g vgq_g];
    inputVec         = [us_vec ug_vec];
end

switch PSS_type
    case 'PSS_on'
        stateVec         = [igd igq phi_fd phi_1d phi_1q phi_2q Dwr theta Pm Dw1 v1 v2 vpss e1 e2 efd e3 ];
    case 'PSS_off'
        stateVec         = [igd igq phi_fd phi_1d phi_1q phi_2q Dwr theta Pm e1 e2 efd e3];
end

% Equilibrium vectors
for i = 1:length(inputVec)
    inputVecEq{i} = sym(strcat(char(inputVec(i)),'_0'));
end
inputVecEq = cell2sym(inputVecEq);

for i = 1:length(stateVec)
    stateVecEq{i} = sym(strcat(char(stateVec(i)),'_0'));
end
stateVecEq = cell2sym(stateVecEq);

for i = 1:length(algVec)
    algVecEq{i} = sym(strcat(char(algVec(i)),'_0'));
end
algVecEq = cell2sym(algVecEq);

%%%%%%%%%%%%%%%%%%%%%%% Differential equations %%%%%%%%%%%%%%%%%%%%%%%%%%%   
% The physical layer (Transformer + Terminal Capacitor):
digd          = (wb/Lt)*(ved - vgd - Rt*igd + wr*Lt*igq);
digq          = (wb/Lt)*(veq - vgq - Rt*igq - wr*Lt*igd);

% SG - Electrical Model 
dphi_fd          = (wb*Rfd/Lad)*efd - (wb*Rfd/Lffd)*(phi_fd - phi_d);
dphi_1d          = (-wb*R1d/L11d)*(phi_1d - phi_d);
dphi_1q          = (-wb*R1q/L11q)*(phi_1q - phi_q);
dphi_2q          = (-wb*R2q/L22q)*(phi_2q - phi_q);

% Mechanical Model: 
dDwr             = (Cm - Ce - KD*Dwr)/(2*H);
% dDelta           = wb*Dwr;
dtheta           = wb*wr;

% Frequency Droop Control
dPm              = (p_ref - DP - Pm)/TG;

% Power System Stabilizer (PSS): 
switch PSS_type
    case 'PSS_on'
        dDw1             = (1/T_LP)*(Dwr - Dw1);
        dv1              = (1/T_HP)*(T_HP*K_PSS*dDw1 - v1);
        dv2              = (1/T1d)*(T1n*dv1 + v1 - v2);
        dvpss            = (1/T2d)*(T2n*dv2 + v2 - vpss);
    case 'PSS_off'
        vpss            = 0;
end

% Excitation Voltage Control:
de1              = (Et - e1)/Tr;
de2              = (Ka*(v_ref - e1 - e3 + vpss) - e2)/Ta;
defd             = (Ke*e2 - efd)/Te;
de3              = (Kfd*defd - e3)/Tfd;

switch PSS_type
    case 'PSS_on'
        diffeqVec = [digd digq dphi_fd dphi_1d dphi_1q dphi_2q dDwr dtheta dPm dDw1 dv1 dv2 dvpss de1 de2 defd de3];    
    case 'PSS_off'
        diffeqVec = [digd digq dphi_fd dphi_1d dphi_1q dphi_2q dDwr dtheta dPm de1 de2 defd de3];    
end

%%%%%%%%%%%%%%%%%%%%%%% Algebraic equations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stator and Rotor currents (Alg.):
Lds = inv(1/Lad + 1/Lffd + 1/L11d);
Lqs = inv(1/Laq + 1/L11q + 1/L22q);

alg1 = wr - Dwr - w_ref;
alg2 = id - ved/Rg - igd;
alg3 = iq - veq/Rg - igq;
alg4 = ved + Ra*id + wr*phi_q;
alg5 = veq + Ra*iq - wr*phi_d;
alg6 = phi_d - Lds*(-id + phi_fd/Lffd + phi_1d/L11d);
alg7 = phi_q - Lqs*(-iq + phi_1q/L11q + phi_2q/L22q);

if i_DG == 1 % If SM is slack (DG=1) no need to convert from global ref. frame to local one (slack is the global one)
    alg8 = vgd - vgd_g; 
    alg9 = vgq - vgq_g;
else         % If SM is NOT slack (DG~=1) we have to convert from global ref. frame to local one (slack is the global one)
    alg8 = vgd - vgd_g*cos(theta_g - theta) + vgq_g*sin(theta_g - theta); 
    alg9 = vgq - vgd_g*sin(theta_g - theta) - vgq_g*cos(theta_g - theta);
end 
alg10 = Cm - Pm/wr;
alg11 = Ce - (phi_d*iq - phi_q*id);
alg12 = DP - (wr-w_ref)/mp;
alg13 = Et - sqrt(ved^2 + veq^2);

algeqVec    = [alg1 alg2 alg3 alg4 alg5 alg6 alg7 alg8 alg9 alg10 alg11 alg12 alg13];

%%%%%%%%%%%%%%%%%%%%%%% Output equations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output vector
if i_DG == 1
    igd_g       = igd;
    igq_g       = igq;
    ved_g       = ved;
    veq_g       = veq;
    outputeqVec_g = [igd_g igq_g theta wr];
else 
    igd_g       = igd*cos(theta - theta_g) - igq*sin(theta - theta_g);
    igq_g       = igd*sin(theta - theta_g) + igq*cos(theta - theta_g);
    ved_g       = ved*cos(theta - theta_g) - veq*sin(theta - theta_g);
    veq_g       = ved*sin(theta - theta_g) + veq*cos(theta - theta_g);    
    outputeqVec_g = [igd_g igq_g]; 
end
p_e              = (ved*igd + veq*igq);
q_e              = (-ved*igq + veq*igd);
outputeqVec_s    = [p_e q_e wr theta dDwr Et];
outputeqVec      = [outputeqVec_s outputeqVec_g];

%%%%%%%%%%%%%%%%%%%%%%% Linearization %%%%%%%%%%%%%%%%%%%%%%%%%%%

Fx          = simplify(jacobian(diffeqVec, stateVec));
Fz          = simplify(jacobian(diffeqVec, algVec));
Fu          = simplify(jacobian(diffeqVec, inputVec));

Gx          = simplify(jacobian(algeqVec, stateVec));
Gz          = simplify(jacobian(algeqVec, algVec));
Gu          = simplify(jacobian(algeqVec, inputVec));
invGz       = inv(Gz);

Hx          = simplify(jacobian(outputeqVec, stateVec));
Hz          = simplify(jacobian(outputeqVec, algVec));
Hu          = simplify(jacobian(outputeqVec, inputVec));

Ai          = Fx - Fz*inv(Gz)*Gx;
Bi          = Fu - Fz*inv(Gz)*Gu;
Ci          = Hx - Hz*inv(Gz)*Gx;
Di          = Hu - Hz*inv(Gz)*Gu;

A_SG0      = subs(Ai, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
B_SG0      = subs(Bi, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
C_SG0      = subs(Ci, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
D_SG0      = subs(Di, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);

switch PSS_type
    case 'PSS_on'
        states_SG  = ["i_{gd}" "i_{gq}" "\Psi_{fd}" "\Psi_{1d}" "\Psi_{1q}" "\Psi_{2q}" ...
                      "\Delta \omega_r" "\theta" "P_m" "\Delta\omega_1" "v_1" "v_2" "v_{PSS}" ...
                      "e_1" "e_2" "e_{fd}" "e_3"];
    case 'PSS_off'
        states_SG  = ["i_{gd}" "i_{gq}" "\Psi_{fd}" "\Psi_{1d}" "\Psi_{1q}" "\Psi_{2q}" ...
                      "\Delta \omega_r" "\theta" "P_m" "e_1" "e_2" "e_{fd}" "e_3"];
end

inputs_SG   = ["V^*" "P^*" "\omega^*"]; % The local inputs names 
outputs_SG   = ["p_e" "q_e" "\omega_r" "\theta" "\Delta\dot\omega_r" "V_t"]; % The local outputs names

nStates_SG  = size(A_SG0,1); 
nIP_SG      = [size(B_SG0,2) length(us_vec) length(ug_vec)]; 
nOP_SG      = [size(C_SG0,1) length(outputeqVec_s) length(outputeqVec_g)]; 

save(strcat('C:\_UsersDatas\keladaf\Library_V4_2024','\Symbolic\ssSG_matrices.mat'),'A_SG0','B_SG0','C_SG0','D_SG0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_SG0./A_SG0); % is zero if (i,j) = true
fid = fopen(strcat('C:\_UsersDatas\keladaf\Library_V4_2024','\Symbolic\A_SG0.txt'),'wt');
for i = 1:length(A_SG0)
    for j = 1:length(A_SG0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_SG0(i,j)));            
        end
    end
end
fclose(fid);

end