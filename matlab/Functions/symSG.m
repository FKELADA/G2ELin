function [A_SG0,B_SG0,C_SG0,D_SG0,states_SG,inputs_SG,outputs_SG,nStates_SG,nIP_SG,nOP_SG] = symSG(PSS,i_DG)

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
syms wb Rt Lt Rg Ra Ll Lad Laq Lfd Rfd L1d L1q L2q R1d R1q R2q H KD mp TG T_LP T_HP K_PSS T1n T1d T2n T2d Tr Ka Ta Ke Te Kfd Tfd

% States (x)
syms igd igq phi_d phi_q phi_fd phi_1d phi_1q phi_2q Dwr theta Pm e1 e2 efd e3 Dw1 v1 v2 vpss

%Algebraic (z)
syms wr ved veq id i1d ifd iq i1q i2q vgd vgq Cm Ce DP Et

% Inputs (u)
syms v_ref p_ref w_ref theta_g vgd_g vgq_g                  %until vdc_ref --> u_k^s (closed loop inputs), after: u^g

%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

algVec           = [wr ved veq id i1d ifd iq i1q i2q vgd vgq Cm Ce DP Et];
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
        stateVec         = [igd igq phi_d phi_q phi_fd phi_1d phi_1q phi_2q Dwr theta Pm Dw1 v1 v2 vpss e1 e2 efd e3 ];
    case 'PSS_off'
        stateVec         = [igd igq phi_d phi_q phi_fd phi_1d phi_1q phi_2q Dwr theta Pm e1 e2 efd e3];
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
dphi_d           = wb*(ved + Ra*id + wr*phi_q);
dphi_q           = wb*(veq + Ra*iq - wr*phi_d);
dphi_fd          = (wb*Rfd/Lad)*efd - wb*Rfd*ifd;
dphi_1d          = wb*(-R1d*i1d);
dphi_1q          = wb*(-R1q*i1q);
dphi_2q          = wb*(-R2q*i2q);

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
        diffeqVec = [digd digq dphi_d dphi_q dphi_fd dphi_1d dphi_1q dphi_2q dDwr dtheta dPm dDw1 dv1 dv2 dvpss de1 de2 defd de3];    
    case 'PSS_off'
        diffeqVec = [digd digq dphi_d dphi_q dphi_fd dphi_1d dphi_1q dphi_2q dDwr dtheta dPm de1 de2 defd de3];    
end

%%%%%%%%%%%%%%%%%%%%%%% Algebraic equations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stator and Rotor currents (Alg.):
invd             = inv([-Lad-Ll Lad Lad; -Lad L1d+Lad Lad; -Lad Lad Lfd+Lad]);
invq             = inv([-Laq-Ll Laq Laq; -Laq L1q+Laq Laq; -Laq Laq L2q+Laq]);
currents_d       = invd * [phi_d phi_1d phi_fd].';
currents_q       = invq * [phi_q phi_1q phi_2q].';

alg1 = wr - Dwr - w_ref;
alg2 = ved - (id - igd)*Rg;
alg3 = veq - (iq - igq)*Rg;
alg4 = id - currents_d(1);
alg5 = i1d - currents_d(2);
alg6 = ifd - currents_d(3);
alg7 = iq - currents_q(1);
alg8 = i1q - currents_q(2);
alg9 = i2q - currents_q(3);
if i_DG == 1 % If SM is slack (DG=1) no need to convert from global ref. frame to local one (slack is the global one)
    alg10 = vgd - vgd_g; 
    alg11 = vgq - vgq_g;
else         % If SM is NOT slack (DG~=1) we have to convert from global ref. frame to local one (slack is the global one)
    alg10 = vgd - vgd_g*cos(theta_g - theta) + vgq_g*sin(theta_g - theta); 
    alg11 = vgq - vgd_g*sin(theta_g - theta) - vgq_g*cos(theta_g - theta);
end 
alg12 = Cm - Pm/wr;
alg13 = Ce - (phi_d*iq - phi_q*id);
alg14 = DP - (wr-w_ref)/mp;
alg15 = Et - sqrt(ved^2 + veq^2);

algeqVec    = [alg1 alg2 alg3 alg4 alg5 alg6 alg7 alg8 alg9 alg10 alg11 alg12 alg13 alg14 alg15];

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
        states_SG  = ["i_{{gd}" "i_{{gq}" "\Psi_{d" "\Psi_{q" "\Psi_{{fd}" "\Psi_{{1d}" "\Psi_{{1q}" "\Psi_{{2q}" ...
                      "\Delta \omega_{r" "\theta_{" "P_{m" "\Delta\omega_{1" "v_{1" "v_{2" "v_{{PSS}" ...
                      "e_{1" "e_{2" "e_{{fd}" "e_{3"];
    case 'PSS_off'
        states_SG  = ["i_{{gd}" "i_{{gq}" "\Psi_{d" "\Psi_{q" "\Psi_{{fd}" "\Psi_{{1d}" "\Psi_{{1q}" "\Psi_{{2q}" ...
                      "\Delta \omega_{r" "\theta_{" "P_{m" "e_{1" "e_{2" "e_{{fd}" "e_{3"];
end

inputs_SG   = ["V^*" "P^*" "\omega^*"]; % The local inputs names 
outputs_SG   = ["p_e" "q_e" "\omega_r" "\theta" "\Delta\dot\omega_r" "V_t"]; % The local outputs names

nStates_SG  = size(A_SG0,1); 
nIP_SG      = [size(B_SG0,2) length(us_vec) length(ug_vec)]; 
nOP_SG      = [size(C_SG0,1) length(outputeqVec_s) length(outputeqVec_g)]; 

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssSG_matrices.mat'),'A_SG0','B_SG0','C_SG0','D_SG0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_SG0./A_SG0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_SG0.txt'),'wt');
for i = 1:length(A_SG0)
    for j = 1:length(A_SG0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_SG0(i,j)));            
        end
    end
end
fclose(fid);

end