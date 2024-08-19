function [A_GFL0,B_GFL0,C_GFL0,D_GFL0,states_GFL,inputs_GFL,outputs_GFL,nStates_GFL,nIP_GFL,nOP_GFL] = symGFL()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GFL #1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% Define symbols %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(strcat(cd,'/Symbolic'));

% Small-signal model for GFL+Transformer 
% DC side considered 

% Parameters
syms wb wff Rf Lf Cf Rt Lt Kpd Kid Kiq KpCL KiCL Kffv Cdc Gdc Tdc Kipll Kppll

% States (x)
syms M_d M_q M_CLd M_CLq M_pll theta_pll isd isq igd igq ved veq vdc idc

% Algebraic (z)
syms md mq w_pll

% Inputs (u)
syms idc_ref vdc_ref q_ref theta_g vgd_g vgq_g                  %until vdc_ref --> u_k^s (closed loop inputs), after: u^g


%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

stateVec         = [isd isq igd igq ved veq vdc idc M_d M_q M_CLd M_CLq M_pll theta_pll];
algVec           = [md mq w_pll];
us_vec           = [vdc_ref q_ref idc_ref];
ug_vec           = [theta_g vgd_g vgq_g];
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

for i = 1:length(algVec)
    algVecEq{i} = sym(strcat(char(algVec(i)),'_0'));
end
algVecEq = cell2sym(algVecEq);
% algVecEq = [];

%%%%%%%%%%%%%%%%% The Physical layer (LCL-Filter + transformer): %%%%%%%%%%%%%%%%%%%%%%

vgd         = vgd_g*cos(theta_g - theta_pll) - vgq_g*sin(theta_g - theta_pll); % have to adjust for GFM-VSC equations
vgq         = vgd_g*sin(theta_g - theta_pll) + vgq_g*cos(theta_g - theta_pll);

disd        = (wb/Lf)*(md*vdc - ved - Rf*isd + w_pll*Lf*isq);
disq        = (wb/Lf)*(mq*vdc - veq - Rf*isq - w_pll*Lf*isd);

digd        = (wb/Lt)*(ved - vgd - Rt*igd + w_pll*Lt*igq);
digq        = (wb/Lt)*(veq - vgq - Rt*igq - w_pll*Lt*igd);

dved        = (wb/Cf)*(isd - igd + w_pll*Cf*veq);
dveq        = (wb/Cf)*(isq - igq - w_pll*Cf*ved);

dvdc        = (wb/Cdc)*(idc - Gdc*vdc - md*isd - mq*isq);

diffeqVec_PL = [disd disq digd digq dved dveq dvdc];

%%%%%%%%%%%%%%%%%%%%%% DC Side %%%%%%%%%%%%%%%%%%%%%%%

didc        = (1/Tdc)*(idc_ref - idc);
diffeqVec_dc = [didc];

%%%%%%%%%%%%%%%% The Outer dq-axis control: %%%%%%%%%%%%%%%

% vdc_ref control (d-axis)
dM_d     = Kid*(vdc_ref - vdc);
isd_ref  = Kpd*(vdc_ref - vdc) + M_d;

% Q_ref loop: 
q        = -ved*igq + veq*igd;
p        = ved*igd + veq*igq;
dM_q     = Kiq*(q_ref - q);

diffeqVec_external = [dM_d dM_q];

%%%%%%%%%%%%%%%%%% The inner current control loop + PLL: %%%%%%%%%%%%%%%%%%%%%

dM_CLd      = KiCL*(isd_ref - isd);
dM_CLq      = KiCL*(M_q - isq);

%%%%%%%%%%%%%%%%%% The PLL: %%%%%%%%%%%%%%%%%%%%%
dM_pll      = Kipll*veq;
dtheta_pll  = wb*w_pll;

diffeqVec_internal = [dM_CLd dM_CLq dM_pll dtheta_pll];

diffeqVec   = [diffeqVec_PL diffeqVec_dc diffeqVec_external diffeqVec_internal];

%%%%%%%%%%%%%%%%%%%%%%% Algebraic equations %%%%%%%%%%%%%%%%%%%%%%%%%%%

alg1        = md - (1/vdc)*(KpCL*(isd_ref - isd)  + Kffv*ved - wff*Lf*isq + M_CLd);
alg2        = mq - (1/vdc)*(KpCL*(M_q - isq)  + Kffv*veq + wff*Lf*isd + M_CLq);
alg3        = w_pll - M_pll - Kppll*veq - wff;
algeqVec    = [alg1 alg2 alg3];

%%%%%%%%%%%%%%%%%%%%%%% Output vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%

igd_g            = igd*cos(theta_pll - theta_g) - igq*sin(theta_pll - theta_g);
igq_g            = igd*sin(theta_pll - theta_g) + igq*cos(theta_pll - theta_g);
outputeqVec_g    = [igd_g igq_g];

w_pll            = M_pll + Kppll*veq + wff;
Vt = sqrt(ved^2 + veq^2);
outputeqVec_s    = [p q w_pll Vt];
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

A_GFL0      = subs(Ai, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
B_GFL0      = subs(Bi, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
C_GFL0      = subs(Ci, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
D_GFL0      = subs(Di, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);

states_GFL  = ["i_{{sd}" "i_{{sq}" "i_{{gd}" "i_{{gq}" "v_{{ed}" "v_{{eq}" "v_{{dc}" "i_{{dc}" ...
               "M_{d" "M_{q" "M_{{{CL}_d}" "M_{{{CL}_q}" "M_{{pll}" "\theta_{{pll}"];
inputs_GFL   = ["V^*_{dc}" "Q^*" "I^*_{dc}"]; % The local inputs names
outputs_GFL   = ["p_e" "q_e" "\omega_{pll}" "V_t"]; % The local outputs names
nStates_GFL = size(A_GFL0,1); 
nIP_GFL     = [size(B_GFL0,2) length(us_vec) length(ug_vec)]; 
nOP_GFL     = [size(C_GFL0,1) length(outputeqVec_s) length(outputeqVec_g)]; 

global baseDirectory
save(strcat(baseDirectory,'\Symbolic\ssGFL_matrices.mat'),'A_GFL0','B_GFL0','C_GFL0','D_GFL0');

%%%%%%%%%%%%%%%%%%%%%%%%%% Print to text file %%%%%%%%%%%%%%%%%%%%%%%%%%%

checkZero = isnan(A_GFL0./A_GFL0); % is zero if (i,j) = true
fid = fopen(strcat(baseDirectory,'\Symbolic\A_GFL0.txt'),'wt');
for i = 1:length(A_GFL0)
    for j = 1:length(A_GFL0)
        if (checkZero(i,j) == false)            
            fprintf(fid, 'A[%d][%d] = %s \n', i, j, char(A_GFL0(i,j)));            
        end
    end
end
fclose(fid);

end