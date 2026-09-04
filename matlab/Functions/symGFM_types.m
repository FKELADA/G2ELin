function [A_GFM0,B_GFM0,C_GFM0,D_GFM0,states_GFM,inputs_GFM,outputs_GFM,nStates_GFM,nIP_GFM,nOP_GFM] = symGFM_types(P_control,i_DG)
% clc 
% clear all 
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GFM #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(strcat(cd,'/Symbolic'));

% Outer power control type:
P_type = P_control;

% Small-signal model for GFM+Transformer 
% DC side considered 

% Parameters
syms wb wff Rf Lf Cf Rt Lt mp nq wf wc KpVL KiVL Kffi KpCL KiCL Kffv Cdc Gdc Kpdc Tdc eta alfa J Dp K Dq K_theta KiPLL KpPLL

% States (x)
syms pm Dw theta qm M_VLd M_VLq M_CLd M_CLq isd isq igd igq ved veq vdc idc ved_ref Phi vdc_m M_PLL 

% Algebraic (z)
syms md mq w

% Inputs (u)
syms p_ref q_ref ve_ref w_ref vdc_ref theta_g vgd_g vgq_g                   %until vdc_ref --> u_k^s (closed loop inputs), after: u^g

%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%
algVec           = [md mq w];
us_vec           = [p_ref q_ref ve_ref w_ref vdc_ref];

if i_DG == 1 
    ug_vec           = [vgd_g vgq_g];
    inputVec         = [us_vec ug_vec];
else 
    ug_vec           = [theta_g vgd_g vgq_g];
    inputVec         = [us_vec ug_vec];
end

switch P_type
    case 'Droop'
        stateVec         = [isd isq igd igq ved veq vdc idc pm theta qm M_VLd M_VLq M_CLd M_CLq];
    case 'Droop+filter'
        stateVec         = [isd isq igd igq ved veq vdc idc pm Dw theta qm M_VLd M_VLq M_CLd M_CLq];
    case 'Droop+PLL'
        stateVec         = [isd isq igd igq ved veq vdc idc pm M_PLL theta qm M_VLd M_VLq M_CLd M_CLq];
    case 'dVOC'
        stateVec         = [isd isq igd igq ved veq vdc idc pm theta qm ved_ref M_VLd M_VLq M_CLd M_CLq];
    case 'VSM'
        stateVec         = [isd isq igd igq ved veq vdc idc Dw theta Phi M_VLd M_VLq M_CLd M_CLq];
    case 'Matching'
        stateVec         = [isd isq igd igq ved veq vdc idc vdc_m theta M_VLd M_VLq M_CLd M_CLq];
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
% algVecEq = [];

%%%%%%%%%%%%%%%%% The Physical layer (LCL-Filter + transformer): %%%%%%%%%%%%%%%%%%%%%%

if i_DG == 1
    vgd         = vgd_g; % have to adjust for GFM-VSC equations
    vgq         = vgq_g;
else
    vgd         = vgd_g*cos(theta_g - theta) - vgq_g*sin(theta_g - theta); % have to adjust for GFM-VSC equations
    vgq         = vgd_g*sin(theta_g - theta) + vgq_g*cos(theta_g - theta);
end    

disd        = (wb/Lf)*(md*vdc - ved - Rf*isd + w*Lf*isq);
disq        = (wb/Lf)*(mq*vdc - veq - Rf*isq - w*Lf*isd);

digd        = (wb/Lt)*(ved - vgd - Rt*igd + w*Lt*igq);
digq        = (wb/Lt)*(veq - vgq - Rt*igq - w*Lt*igd);

dved        = (wb/Cf)*(isd - igd + w*Cf*veq);
dveq        = (wb/Cf)*(isq - igq - w*Cf*ved);

dvdc        = (wb/Cdc)*(idc - Gdc*vdc - md*isd - mq*isq);

diffeqVec_PL = [disd disq digd digq dved dveq dvdc];

%%%%%%%%%%%%%%%%%%%%%% DC Voltage Control %%%%%%%%%%%%%%%%%%%%%%%

idc_ref     = (p_ref/vdc_ref) + Kpdc*(vdc_ref - vdc);
didc        = (idc_ref - idc)/Tdc;

diffeqVec_dc = [didc];

%%%%%%%%%%%%%%%% The Outer power control loop - GFM-types: %%%%%%%%%%%%%%% 
q        = -ved*igq + veq*igd;
p        = ved*igd + veq*igq;

switch P_type
    case 'Droop'
        Dw       = mp*(p_ref - pm);
        dpm      = wf*(p - pm);
        dqm      = wf*(q - qm);
        ved_ref  = ve_ref + (q_ref - qm)*nq;
        veq_ref  = 0;
        dtheta   = wb*w;
        diffeqVec_external = [dpm dtheta dqm];
    case 'Droop+filter'
        dpm      = wf*(p - pm);
        dDw      = mp*wc*(p_ref - pm) - wc*Dw;
        dqm      = wf*(q - qm); 
        ved_ref  = ve_ref + (q_ref - qm)*nq;
        veq_ref  = 0;        
        dtheta   = wb*w;
        diffeqVec_external = [dpm dDw dtheta dqm];
    case 'Droop+PLL'
        Dw       = mp*(p_ref - pm);
        dpm      = wf*(p - pm);
        dqm      = wf*(q - qm);
        dM_PLL   = veq*KiPLL;
        ved_ref  = ve_ref + (q_ref - qm)*nq;
        veq_ref  = 0;
        dtheta   = wb*w;
        diffeqVec_external = [dpm dM_PLL dtheta dqm];
    case 'dVOC'
        Dw       = eta*((p_ref/ve_ref^2) - (pm/ved_ref^2));
        dpm      = wf*(p - pm);
        dqm      = wf*(q - qm);
        dved_ref = (eta*((q_ref/ve_ref^2) - (qm/ved_ref^2)) + (eta*alfa/ve_ref^2)*(ve_ref^2 - ved_ref^2))*ved_ref;
        veq_ref  = 0;
        dtheta   = wb*w;
        diffeqVec_external = [dpm dtheta dqm dved_ref];    
    case 'VSM'
        dDw      = (1/(J*wff))*(p_ref - p) - (Dp/J)*Dw;
        dPhi     = (1/K)*(q_ref - q) + (Dq/K)*(ve_ref - ved);
%         va       = Phi*w*sin(theta);
%         vb       = Phi*w*sin(theta - 2*pi/3);
%         vc       = Phi*w*sin(theta + 2*pi/3);
%         ved_ref  = (2/3)*(va*sin(theta) + vb*sin(theta - 2*pi/3) + vc*sin(theta + 2*pi/3));
        ved_ref  = w*Phi;
        veq_ref  = 0;
        dtheta   = wb*w;      
        diffeqVec_external = [dDw dtheta dPhi];  
    case 'Matching'
        Dw       = K_theta*(vdc_m - vdc_ref);
        dvdc_m   = wf*(vdc - vdc_m);
        ved_ref  = ve_ref;
        veq_ref  = 0;
        dtheta   = wb*w;
        diffeqVec_external = [dvdc_m dtheta]; 
end

%%%%%%%%%%%%%%%%%% The inner voltage control loop: %%%%%%%%%%%%%%%%%%%%
      
dM_VLd      = KiVL*(ved_ref - ved);
dM_VLq      = KiVL*(veq_ref - veq);

%%%%%%%%%%%%%%%%%% The inner current control loop: %%%%%%%%%%%%%%%%%%%%%

isd_ref     = KpVL*(ved_ref - ved) + M_VLd + Kffi*igd - wff*Cf*veq;
isq_ref     = KpVL*(veq_ref - veq) + M_VLq + Kffi*igq + wff*Cf*ved;
dM_CLd      = KiCL*(isd_ref - isd);
dM_CLq      = KiCL*(isq_ref - isq);

diffeqVec_internal = [dM_VLd dM_VLq dM_CLd dM_CLq];

%%%%%%%%%%%%%%%%%%%%%%% Linearization %%%%%%%%%%%%%%%%%%%%%%%%%%%
diffeqVec   = [diffeqVec_PL diffeqVec_dc diffeqVec_external diffeqVec_internal];

% Algebraic Equations:
alg1        = md - (1/vdc)*(KpCL*(isd_ref - isd)  + Kffv*ved - wff*Lf*isq + M_CLd);
alg2        = mq - (1/vdc)*(KpCL*(isq_ref - isq)  + Kffv*veq + wff*Lf*isd + M_CLq);
switch P_type
    case 'Droop'
        alg3        = w - Dw - w_ref;
    case 'Droop+filter'
        alg3        = w - Dw - w_ref;
    case 'Droop+PLL'
        w_PLL       = wff + KpPLL*veq + M_PLL;
        alg3        = w - Dw - w_PLL;
    case 'dVOC'
        alg3        = w - Dw - w_ref;
    case 'VSM'
        alg3        = w - Dw - w_ref;
    case 'Matching'
        alg3        = w - Dw - w_ref;
end

algeqVec    = [alg1 alg2 alg3];


% Output vector
if i_DG == 1
    igd_g       = igd;
    igq_g       = igq;
    outputeqVec_g = [igd_g igq_g theta w];
else 
    igd_g       = igd*cos(theta - theta_g) - igq*sin(theta - theta_g);
    igq_g       = igd*sin(theta - theta_g) + igq*cos(theta - theta_g);
    outputeqVec_g = [igd_g igq_g];
end

Vt = sqrt(ved^2 + veq^2);
outputeqVec_s    = [p q w Vt];
outputeqVec      = [outputeqVec_s outputeqVec_g];

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

A_GFM0      = subs(Ai, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
B_GFM0      = subs(Bi, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
C_GFM0      = subs(Ci, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
D_GFM0      = subs(Di, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);

switch P_type
    case 'Droop'
        states_GFM  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
                       "p_{m" "\theta_{" "q_{m" "M_{{{VL}_d}" "M_{{{VL}_q}" "M_{{{CL}_d}" "M_{{{CL}_q}"];
    case 'Droop+filter'
        states_GFM  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
                       "p_{m" "\delta\omega_{" "\theta_{" "q_{m" "M_{{{VL}_d}" "M_{{{VL}_q}" "M_{{{CL}_d}" "M_{{{CL}_q}"];
    case 'Droop+PLL'
        states_GFM  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
                       "p_{m" "M_{pll}" "\theta_{" "q_{m" "M_{{{VL}_d}" "M_{{{VL}_q}" "M_{{{CL}_d}" "M_{{{CL}_q}"];
    case 'dVOC'
        states_GFM  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
                       "p_{m" "\theta_{" "q_{m" "v^*_{{ed}" "M_{{{VL}_d}" "M_{{{VL}_q}" "M_{{{CL}_d}" "M_{{{CL}_q}"];  
    case 'VSM'
        states_GFM  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
                       "\delta\omega_{" "\theta_{" "\phi_{" "M_{{{VL}_d}" "M_{{{VL}_q}" "M_{{{CL}_d}" "M_{{{CL}_q}"]; 
    case 'Matching'
        states_GFM  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
                       "v_{{{dc}_m}" "\theta_{" "M_{{{VL}_d}" "M_{{{VL}_q}" "M_{{{CL}_d}" "M_{{{CL}_q}"]; 
end

inputs_GFM   = ["P^*" "Q^*" "V^*" "\omega^*" "V^*_{dc}"]; % the local input names
outputs_GFM   = ["p_e" "q_e" "\omega" "V_t"]; % the local output names
nStates_GFM = size(A_GFM0,1);
nIP_GFM     = [size(B_GFM0,2) length(us_vec) length(ug_vec)]; 
nOP_GFM     = [size(C_GFM0,1) length(outputeqVec_s) length(outputeqVec_g)];

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssGFM_matrices_',P_type,'.mat'),'A_GFM0','B_GFM0','C_GFM0','D_GFM0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_GFM0./A_GFM0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_GFM0_',P_type,'.txt'),'wt');
for i = 1:length(A_GFM0)
    for j = 1:length(A_GFM0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_GFM0(i,j)));            
        end
    end
end
fclose(fid);

end