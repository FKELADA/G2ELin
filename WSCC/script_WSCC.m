clc; clear ; close all
all_fig = findall(0, 'type', 'figure');
close(all_fig)
warning('off', 'all');

%% ========================================================================
%  Section I: Define simulation/analysis parameters
%  ========================================================================
%% Simulation Parameters (for the user - optional)
simlt.T_start          = 0;        % Simulation start time [s]
simlt.T_stop           = 5;        % Simulation end time [s]
simlt.T_step           = 1e-5;     % Simulation time step (if fixed) [s]
%% User specified inputs:
model.name             = ['WSCC_1SM_1GFM_1GFL'];  % The simulink file name
Dynamic_sim            = 0; %yes=1, no=0 
ss_analysis            = 1; %yes=1, no=0
% Only in case of dynamic_simulations:
reinit_LF              = 0; % In case you want to run Simulink loadflow only: 1: run load flow in simulink for initialisations 
runSim                 = 0; % run the simulation directly when you run this script
% analysis:
Participation          = 1; single_mode = 1; mode_participation   = 3;
EigenVal_Plot          = 1; Balance_option = 'nobalance';
mode_shape             = 1; mode_mshape          = 3;
sensitivity            = 1; mode_sensitivity     = 3;
step_resp              = 0; input     = "P^*_{GFM_1}"; output    = "p_e_s"; StepM = 0.1; %You can check all inputs and outputs in the variables: inputNames and outputNames
Free_motion            = 0; free_resp_state = "\Delta \omega_{r_s}";
%% ========================================================================
%  Section II: Network Description:
%  ========================================================================
% General network parameters:
Parameters.fb           = 60; %in Hz
Parameters.wb           = 2*pi*Parameters.fb; 
Parameters.Sb           = 100e6; % in VA
Parameters.Un_LV        = 18e3; % in V
Parameters.Un_MV        = 60e3; % in V
Parameters.Un_HV        = 230e3; % in V
%----------------------------------------------------------%
% Per unit calculations:
[Ub_LV,Ib_LV,Zb_LV,Lb_LV,Cb_LV,Ub_LV_dc,Ib_LV_dc,Zb_LV_dc,Cb_LV_dc] = PU_calc(Parameters.wb,Parameters.Un_LV,Parameters.Sb);
[Ub_MV,Ib_MV,Zb_MV,Lb_MV,Cb_MV,Ub_MV_dc,Ib_MV_dc,Zb_MV_dc,Cb_MV_dc] = PU_calc(Parameters.wb,Parameters.Un_MV,Parameters.Sb);
[Ub_HV,Ib_HV,Zb_HV,Lb_HV,Cb_HV,Ub_HV_dc,Ib_HV_dc,Zb_HV_dc,Cb_HV_dc] = PU_calc(Parameters.wb,Parameters.Un_HV,Parameters.Sb);
%----------------------------------------------------------%
% Presets:
switch model.name
    case 'WSCC_3SM'
        % 3SM:
        Y_DER     = [1         4        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         4        0     2        1       0            1.05    0       0.001     0       0          0.0625      0;
                     3         4        0     2        1       0            1.05    0       0.001     0       0          0.0586      0];
        
    case 'WSCC_2SM_1GFL'
        % 2SM_1GFL:
        Y_DER     = [1         4        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         4        0     2        1       0            1.05    0       0.001     0       0          0.0625      0;
                     3         3        0     3        1       0            1.05    0       0         0       0          0.0586      0];
        
    case 'WSCC_1SM_2GFL'
        % 1SM_2GFL:
        Y_DER     = [1         4        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         3        0     3        1       0            1.05    0       0         0       0          0.0625      0;
                     3         3        0     3        1       0            1.05    0       0         0       0          0.0586      0];
    case 'WSCC_1SM_1GFM_1GFL'
        % 1SM_1GFM_1GFL         
        Y_DER     = [1         4        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         2        0     2        1       0            1.05    0       0         0       0          0.0625      0;
                     3         3        0     3        1       0            1.05    0       0         0       0          0.0586      0];

    case 'WSCC_1SM_2GFM'
        % 1SM_2GFM         
        Y_DER     = [1         4        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         2        0     2        1       0            1.05    0       0         0       0          0.0625      0;
                     3         2        0     2        1       0            1.05    0       0         0       0          0.0586      0];
    
    case 'WSCC_2SM_1GFM'
        % 2SM_1GFM         
        Y_DER     = [1         4        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     3         2        0     2        1       0            1.05    0       0         0       0          0.0586      0;                     
                     2         4        0     2        1       0            1.05    0       0.001     0       0          0.0625      0.2];
                 
    case 'WSCC_1GFM_2GFL'
        % 1GFM_2GFL        
        Y_DER     = [1         2        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         3        0     3        1       0            1.05    0       0         0       0          0.0625      0;
                     3         3        0     3        1       0            1.05    0       0         0       0          0.0586      0];

    case 'WSCC_3GFM'
        % 3GFM         
        Y_DER     = [1         2        0     1        1.025   0            0.8     0       0.001     0       0          0.0576      0.2;
                     2         2        0     2        1       0            1.05    0       0         0       0          0.0625      0;
                     3         2        0     2        1       0            1.05    0       0         0       0          0.0586      0];

end
% Selection of Different Elements Voltage Level & Nominal Power:
Nom_Power     = [Parameters.Sb; Parameters.Sb; Parameters.Sb; Parameters.Sb; Parameters.Sb; Parameters.Sb]; %1: IB (if any), 2: SMs, 3: TRs, 4: Lines, 5: Loads, 6:VSCs
Voltage_level = ["LV"; "LV"; "HV"; "LV"; "HV"; "HV"; "LV"]; %1: IB (if any), 2: SMs, 3: TRs (prim. side), 4: TRs (sec. side), 5: Lines, 6: Loads, 7: VSCs
[IB_BP,SM_BP,TR1_BP,TR2_BP,Lines_BP,Loads_BP,VSC_BP] = Base_parameters(Voltage_level,Parameters.wb,Parameters.Sb,Parameters.Un_LV,Parameters.Un_MV,Parameters.Un_HV);

% Selection of other specific models:
GFM_P_control = 'Droop'; % GFM controllers available (also changes models in time domain simulations): 'Droop', 'Droop+filter', 'dVOC', 'VSM', 'Matching'
SM_PSS = 'PSS_on'; %PSS on/off (also changes models in time domain simulations): 'PSS_on', 'PSS_off'

% The function "WSCC_raw: describes the WSCC network basic topology without
% any Generation Unit:
[Y_network1,Y_line1,Y_TR1,nodes_from_r,nodes_to_r,Line_RX,Load.P,Load.Q,Line.R,Line.L,Line.C,Line.length] = WSCC_raw(Parameters.wb,Nom_Power(4),Lines_BP(4),Lines_BP(5),1);

% The function "network_form" forms the network+DER positions by combining
% the raw topology structure to the required positions of the DERs:
[Y_network,Y_line,Y_TR,Y_DG] = network_form(Y_DER,Y_network1,Y_line1,Y_TR1,nodes_from_r,nodes_to_r,GFM_P_control);


%% ========================================================================
%  Section III: Load Flow
%  ========================================================================
[n_DG,n_IB,n_GFM,n_GFL,n_SG,n_loads,n_nodes,n_lines,n_TR,sflag_IB,sflag_GFM,sflag_SG] = counting(Y_network,Y_line,Y_TR);
[N_s, N_G, N_CF, N_Cf, N_L] = nodes_set(Y_network,Y_DER);

%Load Flow:
Sb_MW = Parameters.Sb / 1e6;
[bus_sol,P_loss,I_line,iter,conv_flag,Y_bus] = Load_flow(Y_network,Y_line,Sb_MW);

G_graph = draw_net(Y_network, Y_line);
%% ========================================================================
%  Section III: Events Tested
%  ========================================================================
% Tests of SMs:
T_Pref    = [200 200 200];
delta_Pref= 0.1 .* ones(1,n_SG);

T_Vref    = [200 200 200];
delta_Vref= 0.01 .* ones(1,n_SG);

% Tests of Loads:
T_Pload    = 250 .* ones(1,n_nodes);
delta_Pload= 0.05*315e6 .* ones(1,n_nodes);

T_Qload    = 500 .* ones(1,n_nodes);
delta_Qload= 0.01*115e6 .* ones(1,n_nodes);

% Test of islanding:
Tsw_breaker = [80 120];

% Tests of GFMs:
Tsw_Pref_GFM    = 200 .* ones(1,n_GFM);
Delta_Pref_GFM  = -0.1 .* ones(1,n_GFM);

Tsw_Vref_GFM    = 300 .* ones(1,n_GFM);
Delta_Vref_GFM  = -0.01 .* ones(1,n_GFM);

% Tests of GFLs:
Tsw_Pref_GFL    = 100 .* ones(1,n_GFL);
Delta_Pref_GFL  = 0.1 .* ones(1,n_GFL);

Tsw_Qref_GFL    = 200 .* ones(1,n_GFL);
Delta_Qref_GFL  = 0.1 .* ones(1,n_GFL);

%% ========================================================================
%  Section V: Defining Default Tuning for all library models
%  ========================================================================
if n_SG > 0
% Synchronous Machine Parameters:
% Nominal Values:
SM.Sn     = Nom_Power(2);
SM.wr     = 1;
% Base Values:
switch Voltage_level(2)
    case 'HV'
        SM.Un     = Parameters.Un_HV;
        SM.Es_base= Ub_HV;
        SM.Is_base= Ib_HV;
        SM.Cb     = Cb_HV;
        SM.Zb     = Zb_HV;
    case 'MV'
        SM.Un     = Parameters.Un_MV;
        SM.Es_base= Ub_MV;
        SM.Is_base= Ib_MV;
        SM.Cb     = Cb_MV;
        SM.Zb     = Zb_MV;
    case 'LV'
        SM.Un     = Parameters.Un_LV;
        SM.Es_base= Ub_LV;
        SM.Is_base= Ib_LV;
        SM.Cb     = Cb_LV;
        SM.Zb     = Zb_LV;
end
% Electrical Parameters (in per unit):
SM.Ll     = 0.15;                   %Stator leakage inductance
SM.Ra     = 0.003;                  %Stator resistance
SM.Lad    = 1.66;                   %Mutual inductance - d-axis
SM.Laq    = 1.61;                   %Mutual inductance - q-axis
SM.Ld     = SM.Ll + SM.Lad;         %Self inductance - d-axis 
SM.Lq     = SM.Ll + SM.Laq;         %Self inductance - q-axis
SM.Lafd   = 1.66;                   %Field Mutual inductance
SM.Lfd    = 0.165;                  %Field leakage inductance
SM.Lffd   = SM.Lafd + SM.Lfd;       %Field Self inductance
SM.Rfd    = 0.0006;                 %Field resistnace
SM.Lf1d   = SM.Lffd - SM.Lfd;       %Mutual inductance field and amortisseur 1 circuit in d-axis (see Kundur equation 3.135 page 89)
SM.L1d    = 0.1713;                 %Mutual inductance stator and amortisseur 1 circuit in d-axis
SM.R1d    = 0.0284;                 %Resistance of amortisseur 1 circuit - d-axis
SM.L1q    = 0.7252;                 %Self inductance of amortisseur 1 circuit - d-axis
SM.R1q    = 0.00619;                %Resistance of amortisseur 1 circuit - q-axis
SM.L2q    = 0.125;                  %Mutual inductance stator and amortisseur 2 circuit in q-axis
SM.R2q    = 0.02368;                %Resistance of amortisseur 2 circuit - q-axis
SM.L11d   = SM.L1d + SM.Lf1d;       %Self inductance amortisseur 1 circuit in d-axis(see kundur equations 3.136 - 3.138, P. 89)
SM.L11q   = SM.L1q + SM.Laq;        %Self inductance amortisseur 1 circuit in q-axis
SM.L22q   = SM.L2q + SM.Laq;        %Self inductance amortisseur 2 circuit in q-axis
SM.Cg     = 1e-9;                   %USELESS - TO BE REMOVED
SM.Cg_pu  = SM.Cg / SM.Cb;          %USELESS - TO BE REMOVED
SM.PL     = 100e3;                  %load in parallel with the SM (Simulink needs a small parallel load with the SM)
SM.RL     = SM.Un*SM.Un/SM.PL;      %equivalent resistnace of the load in parallel with the SM
SM.RL_pu  = SM.RL/SM.Zb;            %the equivalent resistnace in pu

% Mechanical Parameters:
SM.H      = 5 .* ones(1,n_SG);      %Inertia constant (s)
SM.K_D    = 0 .* ones(1,n_SG);      %Damping coefficient (pu)

% Controls Parameters:
% 1) The TG - (Governor + Droop):
SM.mp     = 0.5 .* ones(1,n_SG);    % Governor droop value in percent (%)
SM.Tg     = 0.2;                    % Turbine time constant in seconds

% 2) The AVR:
SM.Tr     = 20e-3;                  % Low-pass filter time constant Tr(s)
SM.Ka     = 300;                    % Regulator gain and time constant [ Ka()  Ta(s) ]
SM.Ta     = 0.001;
SM.Ke     = 1;                      % Exciter  [ Ke()  Te(s) ]
SM.Te     = 0.0001;
SM.Tb     = 0;                      % Transient gain reduction [ Tb(s)  Tc(s) ]
SM.Tc     = 0;
SM.Kf     = 0.001;                  % Damping filter gain and time constant  [ Kf()  Tf(s)  ]
SM.Tf     = 0.1;
SM.Efmin  = -11.5;                  % Regulator output limits and gain [ Efmin, Efmax (pu), Kp() ]
SM.Efmax  = 11.5; 
SM.Kp     = 0;

% 3) The PSS:
SM.T_LP   = 0.03;                   %Sensor time constant (s)         
SM.K_PSS  = 2;                      %PSS gain 
SM.T_HP   = 2;                      %Wash-out time constant (s) 
SM.T1n    = 0.05;                   %Lead-lag (1) time constants (s) 
SM.T1d    = 0.02;                   %Lead-lag (1) time constants (s) 
SM.T2n    = 3;                      %Lead-lag (2) time constants (s) 
SM.T2d    = 5.4;                    %Lead-lag (2) time constants (s)
SM.VSmin  = -0.15;                  %saturation minimum 
SM.VSmax  = 0.15;                   %saturation maximum

% Saturation parameters:
SM.A_sat = 0.031;
SM.B_sat = 6.93;
SM.Phi_T1 = 0.8;
SM.Phi_T2 = 1.2;

end
%%
%____________________________________________________________%

if n_GFM > 0
% VSC - GFM Parameters:
% Nominal Values:
GFM.Sn    = Nom_Power(6);
% Base Values:
switch Voltage_level(7)
    case 'HV'
        GFM.Un     = Parameters.Un_HV;
        GFM.Ub     = Ub_HV;
        GFM.Ib     = Ib_HV;
        GFM.Zb     = Zb_HV;
        GFM.Lb     = Lb_HV;
        GFM.Cb     = Cb_HV;
        GFM.Zb_dc  = Zb_HV_dc;
        GFM.Ub_dc  = Ub_HV_dc;
        GFM.Ib_dc  = Ib_HV_dc;
        GFM.Cb_dc  = Cb_HV_dc;
    case 'MV'
        GFM.Un     = Parameters.Un_MV;
        GFM.Ub     = Ub_MV;
        GFM.Ib     = Ib_MV;
        GFM.Zb     = Zb_MV;
        GFM.Lb     = Lb_MV;
        GFM.Cb     = Cb_MV;        
        GFM.Ub_dc  = Ub_MV_dc;
        GFM.Ib_dc  = Ib_MV_dc;
        GFM.Zb_dc  = Zb_MV_dc;
        GFM.Cb_dc  = Cb_MV_dc;
    case 'LV'
        GFM.Un     = Parameters.Un_LV;
        GFM.Ub     = Ub_LV;
        GFM.Ib     = Ib_LV;
        GFM.Zb     = Zb_LV;
        GFM.Lb     = Lb_LV;
        GFM.Cb     = Cb_LV;        
        GFM.Ub_dc  = Ub_LV_dc;
        GFM.Ib_dc  = Ib_LV_dc;
        GFM.Zb_dc  = Zb_LV_dc;
        GFM.Cb_dc  = Cb_LV_dc;
end
GFM.Un_dc  = GFM.Ub_dc;
GFM.Vdcref = GFM.Un_dc;


GFM.R1_pu           = 0.01 ;        % 1st Filter Resistance (pu)
GFM.L1_pu           = 0.1 ;         % 1st Filter Inductance (pu)
GFM.R2_pu           = 0.0533 ;      % TO BE REMOVED - USELESS % 2nd Filter Resistance (pu)
GFM.L2_pu           = 6.1115e-4;    % TO BE REMOVED - USELESS % 2nd Filter Inductance (pu)
GFM.C_pu            = 0.1;          % Filter capacitance (pu)

% Enables/disables:
SFR.Enable              = 0 .* ones(1,n_GFM);         %Secondary FRequency control by the VSCs enable=1, disable=0
Inner_current_loop.Kffv = 1 .* ones(1,n_GFM);         %Voltage feed-forward in the current loop enable=1, disable=0
Inner_voltage_loop.Kffi = 1 .* ones(1,n_GFM);         %Current feed-forward in the voltage loop enable=1, disable=0
VI.Enable               = 0 .* ones(1,n_GFM);         %Virtual impedance function enable=1, disable=0
Damping_resistor.Enable = 0 .* ones(1,n_GFM);         %Damping resistance function enable=1, disable=0
Vdc_loop.Kffi           = 1 .* ones(1,n_GFM);         %DC-Current feed-forward in the dc-voltage loop enable=1, disable=0

% Default Parameters of the Converter Model and its controls:
%-----------------------------------------------------------------------%
% Physical Model: 
%-----------------------------------------------------------------------%
DC_link.R               = 500e3;                               %DC-link Resistor (ohm)
DC_link.R_pu            = DC_link.R ./ GFM.Zb_dc;              %DC-link Resistor (p.u)
DC_link.G_pu            = 1 ./ DC_link.R_pu;                   %DC-link conductance (p.u)
DC_link.C               = 50e-3;                               %DC-link Capacitance (F)
DC_link.C_pu            = DC_link.C ./ GFM.Cb_dc;              %DC-link Capacitance (p.u)
DC_link.Vdc_init        = GFM.Un_dc;                           %DC-link capacitor initial voltage (V)
DC_link.T_dc            = 1e-3 .* ones(1,n_GFM);               %DC-source physical time delay (s)
Idc_max                 = 1.4*GFM.Ib_dc;                       %Maximum DC-source current

GFM.R1               = GFM.R1_pu .* GFM.Zb;                    % Filter Resistance (ohms)
GFM.L1               = GFM.L1_pu .* GFM.Lb;                    % Filter Inductance (H)
GFM.R2               = GFM.R2_pu .* GFM.Zb;                    % TO BE REMOVED - USELESS Filter Resistance (ohms)
GFM.L2               = GFM.L2_pu .* GFM.Lb;                    % TO BE REMOVED - USELESS Filter Inductance (H)
GFM.C                = GFM.C_pu .* GFM.Cb;                     % Filter Inductance (F)

%-----------------------------------------------------------------------%
% GFM Models
%-----------------------------------------------------------------------%
H.first_order           = 3;                                   %Equivalent inertia constant (1st order droop) (s)
H.second_order          = 6;                                   %Equivalent inertia constant (2nd order droop) (s)
KD_opt                  = 200;                                 %Equivalent damping factor (pu)
Droop.mp                = (1/KD_opt) .* ones(1,n_GFM);         %Active power droop value (pu)
q_filter.wf             = 1./(2.*Droop.mp.*H.first_order);     %1st order filter cut-off frequency (rad/s)
Droop.wc                = (q_filter.wf./(2.*Droop.mp.*q_filter.wf.*H.second_order - 1)) .* ones(1,n_GFM); %2nd order droop filter cut-off frequency (rad/s)

%Q-V Droop: 
Droop.nq                = 0.0001 .* ones(1,n_GFM);             %Reactive power droop (pu)  

%P-f Droop: 
p_filter.T1             = 0.0111;                              % TO BE REMOVED - USELESS
p_filter.T2             = 0.0333;                              % TO BE REMOVED - USELESS
Droop.Int_ini           = 0;
%-----------------------------------------------------------------------%
% dVOC Controller:
dVOC.eta                = Droop.mp; 
dVOC.alpha              = 1./(2.*Droop.nq); 
%-----------------------------------------------------------------------%
% VSM Model: 
VSM.D                   = 1./Droop.mp;                                                          
VSM.J                   = 1./(Droop.mp.*q_filter.wf); 
VSM.Dq                  = 1./Droop.nq;    
VSM.K                   = 1./(Droop.nq.*q_filter.wf); 
%-----------------------------------------------------------------------%
%DC-Voltage Controller (P-controller in case of Ideal DC-Source) IN PU:
Vdc_loop.k_dc           = 1./DC_link.G_pu;
Vdc_loop.tau_dc         = DC_link.C_pu./(Parameters.wb*DC_link.G_pu);
Vdc_loop.tr_dc          = 5e-3;                                                  %Control loop 5% settling time tuning value (s)                                                                 
Vdc_loop.zita           = 0.707;                                                 %Control loop tuning damping ratio.
Vdc_loop.wn             = 3 ./ (Vdc_loop.zita*Vdc_loop.tr_dc);                   
Vdc_loop.Kp             = ((3*Vdc_loop.tau_dc/Vdc_loop.tr_dc)-1)./Vdc_loop.k_dc; %Loop PI controller P-gain (pu)                         
Vdc_loop.Ki             = 0 ;                                                    %Loop PI controller I-gain (pu)
Vdc_loop.Int_ini        = GFM.Vdcref;                   
%-----------------------------------------------------------------------%
% Matching Model: 
Matching.K_theta        = Droop.mp.*Vdc_loop.Kp;  
%-----------------------------------------------------------------------%
% Secondary Frequency Regulation (SFR):
% (NOT USED ANYMORE - TO BE REMOVED)
SFR.Rc                  = 0.001 .* ones(1,n_GFM);
SFR.dP_max              = 0.3 .* ones(1,n_GFM);

% Tertiary Functions: Damping Resistor: (only used in case of transmission
% networks with high X/R ratios (NOT USED ANYMORE - TO BE REMOVED)
Damping_resistor.zita    = 0.5;                                 
Damping_resistor.Rv      = (Damping_resistor.zita*GFM.L2_pu ./ sqrt(1-Damping_resistor.zita^2)) - GFM.R2_pu;
Damping_resistor.W_RV    = [Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 Parameters.wb/5 ];

% Tertiary Functions: Virtual Impedance: (only in case a different R/X
% ratio is to be emulated - have an effect on voltage S.S. value) 
% (NOT USED ANYMORE - TO BE REMOVED)
VI.rho                   = 10;                   
VI.Lv_pu                 = 0.1 .* ones(1,n_GFM);
VI.Rv_pu                 = VI.Lv_pu ./ VI.rho;

% Tertiary Functions: Transient Virtual Impedance: if triggered limit the
% AC current beyond a certain maximum I_max
% (NOT USED ANYMORE - TO BE REMOVED)
TVI.I_thresh_pu          = 1.2 .* ones(1,n_GFM);
TVI.I_max_pu             = 1.4 .* ones(1,n_GFM);
TVI.Rho_XR               = 5 .* ones(1,n_GFM);
TVI.X_0                  = 0.1 .* ones(1,n_GFM); 
TVI.R_0                  = TVI.X_0 ./ TVI.Rho_XR; 
TVI.K_RV_a               = (TVI.I_max_pu - TVI.I_thresh_pu).^2 .* (1 + (TVI.Rho_XR).^2);
TVI.K_RV_b               = 2*(TVI.I_max_pu - TVI.I_thresh_pu) .* (TVI.R_0 + TVI.Rho_XR.*TVI.X_0);
TVI.K_RV_c               = TVI.R_0.^2 + TVI.X_0.^2 - (1./(TVI.I_max_pu.^2));
TVI.K_RV                 = (-TVI.K_RV_b + sqrt(TVI.K_RV_b.^2 - 4*TVI.K_RV_a.*TVI.K_RV_c))./(2.*TVI.K_RV_a);
%__________________________________________________________________________________

%Inner Voltage Loop:
Inner_voltage_loop.kv    = Parameters.wb./GFM.C_pu;
Inner_voltage_loop.zita  = 0.707;                                                                      %Control loop tuning damping ratio.
Inner_voltage_loop.tr    = 15e-3;                                                                      %Control loop 5% settling time tuning value (s)
Inner_voltage_loop.wn    = 3 ./ (Inner_voltage_loop.tr*Inner_voltage_loop.zita);
Inner_voltage_loop.Kp    = 2*Inner_voltage_loop.zita*Inner_voltage_loop.wn ./ Inner_voltage_loop.kv;   %Loop PI controller P-gain (pu)
Inner_voltage_loop.Ki    = Inner_voltage_loop.wn^2 ./ Inner_voltage_loop.kv;                           %Loop PI controller I-gain (pu) 

% Inner Current Loop:
Inner_current_loop.ki    = 1./GFM.R1_pu;
Inner_current_loop.tau_i = GFM.L1_pu./(Parameters.wb .* GFM.R1_pu);
Inner_current_loop.zita  = 0.707;                                                                      %Control loop tuning damping ratio.
Inner_current_loop.tr    = 0.5e-3;                                                                     %Control loop 5% settling time tuning value (s)
Inner_current_loop.wn    = 3 ./ (Inner_current_loop.tr*Inner_current_loop.zita);
Inner_current_loop.Kp    = (2*Inner_current_loop.zita*Inner_current_loop.wn.*Inner_current_loop.tau_i - 1) ./ Inner_current_loop.ki; %Loop PI controller P-gain (pu)
Inner_current_loop.Ki    = (Inner_current_loop.wn^2 .* Inner_current_loop.tau_i) ./ Inner_current_loop.ki;                           %Loop PI controller I-gain (pu)              
end
%%
if n_GFL > 0
%____________________________________________________________%
% VSC - GFL Parameters:
% Nominal Values:
GFL.Sn    = Nom_Power(6);
% Base Values:
switch Voltage_level(7)
    case 'HV'
        GFL.Un     = Parameters.Un_HV;
        GFL.Ub     = Ub_HV;
        GFL.Ib     = Ib_HV;
        GFL.Zb     = Zb_HV;
        GFL.Lb     = Lb_HV;
        GFL.Cb     = Cb_HV;
        GFL.Zb_dc  = Zb_HV_dc;
        GFL.Ub_dc  = Ub_HV_dc;
        GFL.Ib_dc  = Ib_HV_dc;
        GFL.Cb_dc  = Cb_HV_dc;
    case 'MV'
        GFL.Un     = Parameters.Un_MV;
        GFL.Ub     = Ub_MV;
        GFL.Ib     = Ib_MV;
        GFL.Zb     = Zb_MV;
        GFL.Lb     = Lb_MV;
        GFL.Cb     = Cb_MV;        
        GFL.Ub_dc  = Ub_MV_dc;
        GFL.Ib_dc  = Ib_MV_dc;
        GFL.Zb_dc  = Zb_MV_dc;
        GFL.Cb_dc  = Cb_MV_dc;
    case 'LV'
        GFL.Un     = Parameters.Un_LV;
        GFL.Ub     = Ub_LV;
        GFL.Ib     = Ib_LV;
        GFL.Zb     = Zb_LV;
        GFL.Lb     = Lb_LV;
        GFL.Cb     = Cb_LV;        
        GFL.Ub_dc  = Ub_LV_dc;
        GFL.Ib_dc  = Ib_LV_dc;
        GFL.Zb_dc  = Zb_LV_dc;
        GFL.Cb_dc  = Cb_LV_dc;
end
GFL.Un_dc  = GFL.Ub_dc;
GFL.Vdcref = GFL.Un_dc;


GFL.R1_pu           = 0.015 ;                  % 1st Filter Resistance (pu)
GFL.L1_pu           = 0.1 ;                    % 1st Filter Inductance (pu)
GFL.R2_pu           = 0.0533 ;                 % USELESS - TO BE REMOVED
GFL.L2_pu           = 6.1115e-4;               % USELESS - TO BE REMOVED
GFL.C_pu            = 0.11;                    % Filter capacitance (pu)

% Enables/disables:
GFL_Inner_current_loop.Kffv = 1 .* ones(1,n_GFL);         %Voltage feed-forward in the current loop enable=1, disable=0
GFL_Vdc_loop.Kffi           = 1 .* ones(1,n_GFL);         %DC-Current feed-forward in the dc-voltage loop enable=1, disable=0

% Default Parameters of the Converter Model and its controls:
%-----------------------------------------------------------------------%
% Physical Model: 
%-----------------------------------------------------------------------%
GFL_DC_link.R               = 500e3;                               %DC-link Resistor (ohm)
GFL_DC_link.R_pu            = GFL_DC_link.R ./ GFL.Zb_dc;          %DC-link Resistor (p.u)
GFL_DC_link.G_pu            = 1 ./ GFL_DC_link.R_pu;               %DC-link conductance (p.u)
GFL_DC_link.C               = 50e-3;                               %DC-link Capacitance (F)
GFL_DC_link.C_pu            = GFL_DC_link.C ./ GFL.Cb_dc;          %DC-link Capacitance (p.u)
GFL_DC_link.Vdc_init        = GFL.Un_dc;                           %DC-link capacitor initial voltage (V)
GFL_DC_link.T_dc            = 1e-3 .* ones(1,n_GFL);               %DC-source physical time delay (s)
GFL_Idc_max                 = 1.4*GFL.Ib_dc;                       %Maximum DC-source current

GFL.R1               = GFL.R1_pu .* GFL.Zb;                        % Filter Resistance (ohms)
GFL.L1               = GFL.L1_pu .* GFL.Lb;                        % Filter Inductance (H)
GFL.R2               = GFL.R2_pu .* GFL.Zb;                        % USELESS - TO BE REMOVED
GFL.L2               = GFL.L2_pu .* GFL.Lb;                        % USELESS - TO BE REMOVED
GFL.C                = GFL.C_pu .* GFL.Cb;                         % Filter Inductance (F)

% Outer DQ-Axes Control:
% The d-axis control (Vdc control):
GFL.kd                  = -1./GFL_DC_link.G_pu;
GFL.tau_d               = GFL_DC_link.C_pu./(Parameters.wb .* GFL_DC_link.G_pu);
GFL.zita_d              = 0.707;                                                 %Control loop tuning damping ratio.
GFL.tr_d                = 100e-3;                                                %Control loop 5% settling time tuning value (s)
GFL.wn_d                = 3 ./ (GFL.tr_d*GFL.zita_d);
GFL.vdc_Kp              = (2*GFL.zita_d*GFL.wn_d.*GFL.tau_d - 1) ./ GFL.kd;      %Loop PI controller P-gain (pu)
GFL.vdc_Ki              = (GFL.wn_d^2 .* GFL.tau_d) ./ GFL.kd;                   %Loop PI controller I-gain (pu)

% The q-axis control (Q_ref control):
GFL.tr_q                = 100e-3;                                                %Control loop 5% settling time tuning value (s)
GFL.KQ                  = -3/GFL.tr_q;                                           %Loop Integral controller I-gain (pu)
GFL.KP                  = 20;                                                    %USELESS - TO BE REMOVED

% Inner Current Loop:
GFL_Inner_current_loop.ki    = 1./GFL.R1_pu;
GFL_Inner_current_loop.tau_i = GFL.L1_pu./(Parameters.wb .* GFL.R1_pu);
GFL_Inner_current_loop.zita  = 0.707;                                            %Control loop tuning damping ratio.
GFL_Inner_current_loop.tr    = 10e-3;                                            %Control loop 5% settling time tuning value (s)
GFL_Inner_current_loop.wn    = 3 ./ (GFL_Inner_current_loop.tr*GFL_Inner_current_loop.zita);
GFL_Inner_current_loop.Kp    = (2*GFL_Inner_current_loop.zita*GFL_Inner_current_loop.wn.*GFL_Inner_current_loop.tau_i - 1) ./ GFL_Inner_current_loop.ki; %Loop PI controller P-gain (pu)
GFL_Inner_current_loop.Ki    = (GFL_Inner_current_loop.wn^2 .* GFL_Inner_current_loop.tau_i) ./ GFL_Inner_current_loop.ki;                               %Loop PI controller I-gain (pu)

% PLL:
PLL.k                   = Parameters.wb;
PLL.zita                = 0.707;                                     %Control loop tuning damping ratio.
PLL.tr                  = 50e-3;                                     %Control loop 5% settling time tuning value (s)
PLL.wn                  = 3 ./ (PLL.tr*PLL.zita);
PLL.Kp                  = 2*PLL.zita*PLL.wn ./ PLL.k;                %Loop PI controller P-gain (pu)
PLL.Ki                  = PLL.wn^2 ./ PLL.k;                         %Loop PI controller I-gain (pu)

% VIM: 
VIM.wc                  = 0.05*Parameters.wb;
VIM.wc_ini              = 0;
VIM.Iqd_sat             = 0.25;
VIM.Rr_pu               = 0.0005;
VIM.Lr_pu               = 0.05;
VIM.Lm_pu               = 0.6;
VIM.Rr                  = VIM.Rr_pu * GFL.Zb;
VIM.Lr                  = VIM.Lr_pu * GFL.Lb;
VIM.Lm                  = VIM.Lm_pu * GFL.Lb;
VIM.Kp                  = VIM.Rr ./ VIM.Lr;
VIM.Kd                  = 0.001;
VIM.Nd                  = 100;
VIM.H                   = 5.04;        
VIM.J                   = 2 * VIM.H * GFL.Sn * 1e-6 / (Parameters.wb^2);
VIM.D                   = 0.8;
VIM.Inertia_num         = 1;
VIM.Inertia_den         = [VIM.J , VIM.D];
VIM.Inertia_ini         = 0;
VIM.Ke_num              = 3/2.*VIM.Rr.*VIM.Lm.^2;
VIM.Ke_den              = [VIM.Lr.^2 , VIM.Rr.*VIM.Lr];
VIM.Ke_ini              = 0;
VIM.T_sync              = 0;
VIM.f0                  = 50;
VIM.w0                  = 2*pi*VIM.f0 / Parameters.wb;
end
%%
% Transformer parameters:
if n_TR > 0
TR.Sn  = Nom_Power(3) .* ones(1,n_TR);
switch Voltage_level(3)
    case 'HV'
        TR.Un1 = Parameters.Un_HV;
    case 'MV'
        TR.Un1 = Parameters.Un_MV;      
    case 'LV'
        TR.Un1 = Parameters.Un_LV;
end
switch Voltage_level(4)
    case 'HV'
        TR.Un2 = Parameters.Un_HV;
    case 'MV'
        TR.Un2 = Parameters.Un_MV;      
    case 'LV'
        TR.Un2 = Parameters.Un_LV;
end
TR.R1 = (Y_TR(:,4).')./2; % in pu
TR.L1 = (Y_TR(:,5).')./2; % in pu
TR.R2 = (Y_TR(:,4).')./2; % in pu
TR.L2 = (Y_TR(:,5).')./2; % in pu
TR.Rm = 500; % in pu
TR.Lm = 500; % in pu
end
%%
%____________________________________________________________%

% Loads Parameters:
% Calculation of the equivalent Load_R and Load_X:
if n_loads > 0
switch Voltage_level(6)
    case 'HV'
        Load.Un = Parameters.Un_HV;
        for i = 1:n_nodes-n_DG
            Load.S(i)                 = sqrt((Y_network(i,9)*Parameters.Sb).*(Y_network(i,9)*Parameters.Sb) + (Y_network(i,10)*Parameters.Sb).*(Y_network(i,10)*Parameters.Sb));
            Load.PF(i)                = (Y_network(i,9)*Parameters.Sb) ./ Load.S(i);
            Load.Z(i)                 = (bus_sol(i,2)*Parameters.Un_HV)*(bus_sol(i,2)*Parameters.Un_HV) ./ (Load.S(i));
            Load.theta(i)             = acos(Load.PF(i));
            Load.R(i)                 = Load.Z(i) .* Load.PF(i);
            Load.X(i)                 = Load.Z(i) .* sin(Load.theta(i));
            Load.L(i)                 = Load.X(i) ./ Parameters.wb;
            Load.R_pu(i)              = Load.R(i) ./ Zb_HV;
            Load.L_pu(i)              = Load.L(i) ./ Lb_HV;
        end
    case 'MV'
        Load.Un = Parameters.Un_MV;      
        for i = 1:n_nodes-n_DG
            Load.S(i)                 = sqrt((Y_network(i,9)*Parameters.Sb).*(Y_network(i,9)*Parameters.Sb) + (Y_network(i,10)*Parameters.Sb).*(Y_network(i,10)*Parameters.Sb));
            Load.PF(i)                = (Y_network(i,9)*Parameters.Sb) ./ Load.S(i);
            Load.Z(i)                 = (bus_sol(i,2)*Parameters.Un_MV)*(bus_sol(i,2)*Parameters.Un_MV) ./ (Load.S(i));
            Load.theta(i)             = acos(Load.PF(i));
            Load.R(i)                 = Load.Z(i) .* Load.PF(i);
            Load.X(i)                 = Load.Z(i) .* sin(Load.theta(i));
            Load.L(i)                 = Load.X(i) ./ Parameters.wb;
            Load.R_pu(i)              = Load.R(i) ./ Zb_MV;
            Load.L_pu(i)              = Load.L(i) ./ Lb_MV;
        end
    case 'LV'
        Load.Un = Parameters.Un_LV;
        for i = 1:n_nodes-n_DG
            Load.S(i)                 = sqrt((Y_network(i,9)*Parameters.Sb).*(Y_network(i,9)*Parameters.Sb) + (Y_network(i,10)*Parameters.Sb).*(Y_network(i,10)*Parameters.Sb));
            Load.PF(i)                = (Y_network(i,9)*Parameters.Sb) ./ Load.S(i);
            Load.Z(i)                 = (bus_sol(i,2)*Parameters.Un_LV)*(bus_sol(i,2)*Parameters.Un_LV) ./ (Load.S(i));
            Load.theta(i)             = acos(Load.PF(i));
            Load.R(i)                 = Load.Z(i) .* Load.PF(i);
            Load.X(i)                 = Load.Z(i) .* sin(Load.theta(i));
            Load.L(i)                 = Load.X(i) ./ Parameters.wb;
            Load.R_pu(i)              = Load.R(i) ./ Zb_LV;
            Load.L_pu(i)              = Load.L(i) ./ Lb_LV;
        end
end
Load.R_pu(isnan(Load.R_pu)) = [];
Load.L_pu(isnan(Load.L_pu)) = [];
Load.Sn = Nom_Power(5);
end

%% ========================================================================
%  Section VI: Computing Initializations
%  ========================================================================
% SMs' Initializations:
SM_indices = (find(Y_network(:,2)==4)).';
for i_SG = 1:n_SG
    SM.Pt(i_SG) = bus_sol(SM_indices(i_SG),4);
    SM.Qt(i_SG) = bus_sol(SM_indices(i_SG),5);
    SM.Et(i_SG) = bus_sol(SM_indices(i_SG),2);
    SM.alpha(i_SG) = bus_sol(SM_indices(i_SG),3)*pi/180;
    
    SM.v0(i_SG)      = SM.Et(i_SG)*exp(1i*SM.alpha(i_SG)); %terminal voltage of the SM (Et) in alpha-beta stationary ref. frame
    SM.Ii(i_SG)      = (SM.Pt(i_SG) - 1i*SM.Qt(i_SG))/conj(SM.v0(i_SG));
    SM.It(i_SG)      = abs(SM.Ii(i_SG));
    SM.Phi(i_SG)     = angle(SM.It(i_SG));
    SM.delta0(i_SG)  = angle(SM.v0(i_SG) + (SM.Ra+1i*SM.Lq)*SM.Ii(i_SG));  
    SM.theta_0(i_SG) = SM.delta0(i_SG) - pi/2;
    SM.Ic0(i_SG)     = SM.v0(i_SG)/SM.RL_pu;
    SM.Ig0(i_SG)     = SM.Ii(i_SG) - SM.Ic0(i_SG);
    SM.vdq0(i_SG)    = SM.v0(i_SG)*exp(1i*(SM.alpha(i_SG)-SM.theta_0(i_SG)));
    SM.idq0(i_SG)    = SM.Ii(i_SG)*exp(1i*(SM.Phi(i_SG)-SM.theta_0(i_SG)));
    SM.igdq0(i_SG)   = SM.Ig0(i_SG)*exp(1i*(SM.Phi(i_SG)-SM.theta_0(i_SG)));
    SM.Te0(i_SG)     = (real(SM.vdq0(i_SG))+SM.Ra*real(SM.idq0(i_SG)))*real(SM.idq0(i_SG)) + (imag(SM.vdq0(i_SG))+SM.Ra*imag(SM.idq0(i_SG)))*imag(SM.idq0(i_SG));

    SM.ed_0(i_SG) = real(SM.vdq0(i_SG));
    SM.eq_0(i_SG) = imag(SM.vdq0(i_SG));

    SM.id_0(i_SG) = real(SM.idq0(i_SG));
    SM.iq_0(i_SG) = imag(SM.idq0(i_SG));
    
    SM.igd_0(i_SG) = real(SM.igdq0(i_SG));
    SM.igq_0(i_SG) = imag(SM.igdq0(i_SG));
    
    % Considering saturation:
    SM.Phi_at(i_SG)  = abs(SM.v0(i_SG) + (SM.Ra+1i*SM.Ll)*SM.Ii(i_SG));
    if SM.Phi_at(i_SG) <= SM.Phi_T1
        SM.Phi_I(i_SG) = 0;
    elseif SM.Phi_T1 <= SM.Phi_at(i_SG) <= SM.Phi_T2
        SM.Phi_I(i_SG) = SM.A_sat*exp(SM.B_sat*(SM.Phi_at(i_SG) - SM.Phi_T1));
%     else
%         SM.Phi_I(i_SG) = SM.Phi_G2 + SM.L_ratio*(SM.Phi_at(i_SG) - SM.Phi_T2) - SM.Phi_at(i_SG);
    end
%     SM.Ksd(i_SG) = SM.Phi_at(i_SG)/(SM.Phi_at(i_SG)+SM.Phi_I(i_SG));
%     SM.Ksq(i_SG) = SM.Ksd(i_SG); 
    SM.Ksd(i_SG) = 1; 
    SM.Ksq(i_SG) = SM.Ksd(i_SG); 
    SM.Lads(i_SG) = SM.Lad * SM.Ksd(i_SG); 
    SM.Laqs(i_SG) = SM.Laq * SM.Ksq(i_SG); 
%
    SM.Psi_d0(i_SG) = SM.eq_0(i_SG) + SM.Ra*SM.iq_0(i_SG);
    SM.Psi_q0(i_SG) = -SM.ed_0(i_SG) - SM.Ra*SM.id_0(i_SG);
    SM.ifd_0(i_SG)  = (SM.eq_0(i_SG) + SM.Ra*SM.iq_0(i_SG) + (SM.Lads(i_SG)+SM.Ll)*SM.id_0(i_SG))/SM.Lads(i_SG);
    SM.efd_0(i_SG)  = SM.Lad*SM.ifd_0(i_SG);
    SM.Psi_fd0(i_SG) = (SM.Lads(i_SG) + SM.Lfd)*SM.ifd_0(i_SG) - SM.Lads(i_SG)*SM.id_0(i_SG);
    SM.Psi_1d0(i_SG) = SM.Lads(i_SG) *(SM.ifd_0(i_SG) - SM.id_0(i_SG));
    SM.Psi_1q0(i_SG) = -SM.Laqs(i_SG) *SM.iq_0(i_SG);
    SM.Psi_2q0(i_SG) = -SM.Laqs(i_SG) *SM.iq_0(i_SG);
    SM.Te_0(i_SG)    = SM.Te0(i_SG);
end
% Loads Initializations:
Loads_indices = (find(Y_network(:,3)==1)).';
for i_load = 1:n_loads
    Load.Pn(i_load) = bus_sol(Loads_indices(i_load),6)*Load.Sn; %in W
    Load.Qn(i_load) = bus_sol(Loads_indices(i_load),7)*Load.Sn; %in VAr
    Load.Vo(i_load) = bus_sol(Loads_indices(i_load),2); %in pu
    Load.Vphase(i_load) = bus_sol(Loads_indices(i_load),3); %in deg.
end

% VSC-GFM Initializations:
GFM_indices = (find(Y_network(:,2)==2)).';
for i_GFM = 1:n_GFM
    GFM.Pref(i_GFM) = bus_sol(GFM_indices(i_GFM),4); %in pu
    GFM.Qref(i_GFM) = bus_sol(GFM_indices(i_GFM),5); %in pu
    GFM.Vref(i_GFM) = bus_sol(GFM_indices(i_GFM),2); %in pu
    q_filter.Int_ini(i_GFM) = GFM.Qref(i_GFM);
    p_filter.Int_ini(i_GFM) = GFM.Pref(i_GFM);
    GFM.Theta_Eg0(i_GFM) = bus_sol(GFM_indices(i_GFM),3)*pi/180;
    
    GFM.Eg0(i_GFM) = GFM.Vref(i_GFM)*exp(1i*GFM.Theta_Eg0(i_GFM));
    GFM.Eg0dq(i_GFM) = GFM.Eg0(i_GFM)*exp(-1i*GFM.Theta_Eg0(i_GFM));
    GFM.Ig0(i_GFM) = (GFM.Pref(i_GFM) - 1i*GFM.Qref(i_GFM))/conj(GFM.Eg0(i_GFM));
    GFM.Ig0dq(i_GFM) = GFM.Ig0(i_GFM)*exp(-1i*GFM.Theta_Eg0(i_GFM));
    GFM.Ic0(i_GFM) = GFM.Eg0(i_GFM)*(1i*GFM.C_pu);
    GFM.Ic0dq(i_GFM) = GFM.Ic0(i_GFM)*exp(-1i*GFM.Theta_Eg0(i_GFM));
    GFM.Is0(i_GFM) = GFM.Ig0(i_GFM) + GFM.Ic0(i_GFM);
    GFM.Is0dq(i_GFM) = GFM.Ig0dq(i_GFM) + GFM.Ic0dq(i_GFM);
    GFM.Vm0(i_GFM) = GFM.Eg0(i_GFM) + GFM.Is0(i_GFM)*(GFM.R1_pu+1i*GFM.L1_pu);
    GFM.Vm0dq(i_GFM) = GFM.Vm0(i_GFM)*exp(-1i*GFM.Theta_Eg0(i_GFM));
    GFM.Theta_Vm0(i_GFM) = angle(GFM.Vm0(i_GFM));
%     GFM.Ig0(i_GFM) = (GFM.Pref(i_GFM) - i*GFM.Qref(i_GFM))/conj(GFM.eg0(i_GFM));
%     GFM.It(i_GFM) = abs(GFM.Ig0(i_GFM));
%     GFM.Ic(i_GFM) = GFM.eg0(i_GFM)*(i*GFM.C_pu);
%     GFM.Phi(i_GFM) = angle(GFM.Pref(i_GFM) - i*GFM.Qref(i_GFM));
%     GFM.isdq0(i_GFM) = GFM.Ig0(i_GFM)+GFM.Ic(i_GFM);
%     GFM.vmdq0(i_GFM) = (GFM.eg0(i_GFM)+GFM.isdq0(i_GFM)*((GFM.R1_pu+i*GFM.L1_pu)));
%     GFM.delta0(i_GFM) = angle(GFM.vmdq0(i_GFM));
%     
%     GFM.vmd_0(i_GFM) = real(GFM.vmdq0(i_GFM));
%     GFM.vmq_0(i_GFM) = imag(GFM.vmdq0(i_GFM));
% 
%     GFM.isd_0(i_GFM) = real(GFM.isdq0(i_GFM));
%     GFM.isq_0(i_SG) = imag(GFM.isdq0(i_GFM));

    % Grid voltage
    Eg0_d(i_GFM)  = real(GFM.Eg0dq(i_GFM)); 
    Eg0_q(i_GFM)  = imag(GFM.Eg0dq(i_GFM));
    Eg0_dq(i_GFM,:) = [Eg0_d(i_GFM); Eg0_q(i_GFM); 0]'; %column_1 is Ed, column_2 is Eq and column_3 is E0

    P0(i_GFM) = GFM.Pref(i_GFM);
    Q0(i_GFM) = GFM.Qref(i_GFM);
    Ws0       = Parameters.wb;

    % Ig
    Ig0_d(i_GFM)  = real(GFM.Ig0dq(i_GFM));
    Ig0_q(i_GFM)  = imag(GFM.Ig0dq(i_GFM)); 
    Ig0_dq(i_GFM,:) = [Ig0_d(i_GFM); Ig0_q(i_GFM) ; 0]'; %column_1 is Igd, column_2 is Igq and column_3 is Ig0

    Is0_d(i_GFM)  = real(GFM.Is0dq(i_GFM)); 
    Is0_q(i_GFM)  = imag(GFM.Is0dq(i_GFM)); 
    Is0_dq(i_GFM,:) = [Is0_d(i_GFM); Is0_q(i_GFM); 0]';

    % Vconv
    Vm0_d(i_GFM)  = real(GFM.Vm0dq(i_GFM));
    Vm0_q(i_GFM)  = imag(GFM.Vm0dq(i_GFM));
    Vm0_dq(i_GFM,:)  = [Vm0_d(i_GFM); Vm0_q(i_GFM); 0]';

%     % Vgrid
%     Vg0_d(i_GFM)  = -(GFM.R2_pu .* Ig0_d(i_GFM)) + (GFM.L2_pu * 1 .* Ig0_q(i_GFM)) + Eg0_d(i_GFM);
%     Vg0_q(i_GFM)  = -(GFM.L2_pu .* 1 .* Ig0_d(i_GFM)) - (GFM.R2_pu .* Ig0_q(i_GFM)) + Eg0_q(i_GFM);
%     Vg0_dq(i_GFM,:)  = [Vg0_d(i_GFM); Vg0_q(i_GFM); 0]';

 
    theta0(i_GFM) = GFM.Theta_Eg0(i_GFM);
    M_transform_inv = [   sin(theta0(i_GFM))     ,    cos(theta0(i_GFM))     ,    1     ;
                          sin(theta0(i_GFM)-2*pi/3)     ,  cos(theta0(i_GFM)-2*pi/3)  , 1  ;
                          sin(theta0(i_GFM)+2*pi/3)     ,  cos(theta0(i_GFM)+2*pi/3)  ,  1 ];
%
    Eg0_abc(i_GFM,:) = (M_transform_inv * Eg0_dq(i_GFM,:)')';
    Vm0_abc(i_GFM,:) = (M_transform_inv * Vm0_dq(i_GFM,:)')';
    Ig0_abc(i_GFM,:) = (M_transform_inv * Ig0_dq(i_GFM,:)')';
    Is0_abc(i_GFM,:)= (M_transform_inv * Is0_dq(i_GFM,:)')';


    % Inner Voltage Loop Integrator Initilization: 
    Inner_voltage_loop.Int_ini_d(i_GFM) = Is0_d(i_GFM) - (Ig0_d(i_GFM) .* Inner_voltage_loop.Kffi(i_GFM)) + (Eg0_q(i_GFM) .* GFM.C_pu * 1) ; %M_VLde
    Inner_voltage_loop.Int_ini_q(i_GFM) = Is0_q(i_GFM) - (Ig0_q(i_GFM) .* Inner_voltage_loop.Kffi(i_GFM)) - (Eg0_d(i_GFM) .* GFM.C_pu * 1) ; %M_VLqe

    % Inner Current Loop Integrator Initilization: 
    Inner_current_loop.Int_ini_d(i_GFM) = Vm0_d(i_GFM) - (Eg0_d(i_GFM) .* Inner_current_loop.Kffv(i_GFM)) + (Is0_q(i_GFM) .* GFM.L1_pu * 1) ;
    Inner_current_loop.Int_ini_q(i_GFM) = Vm0_q(i_GFM) - (Eg0_q(i_GFM) .* Inner_current_loop.Kffv(i_GFM)) - (Is0_d(i_GFM) .* GFM.L1_pu * 1) ;

    % Initialization of the Damping Resistor Integrator: 
    % Damping_resistor.Int_ini = Damping_resistor.Rv .* Ig0_dq;
    Damping_resistor.Int_ini(i_GFM) = 0;

    % Filter Currents and voltage initialization:
    Filter.I_ini(i_GFM,:)             = Is0_abc(i_GFM,:).*GFM.Ib ;  %
    Filter.V_ini(i_GFM,:)             = Eg0_abc(i_GFM,:).*GFM.Ub ;  %
    Line.I_ini(i_GFM,:)               = Ig0_abc(i_GFM,:).*GFM.Ib ;  %
    % dVOC initialization:
    dVOC.ini_theta(i_GFM)           = GFM.Theta_Eg0(i_GFM);

    % Matching initialization:
    Matching.Int_ini(i_GFM)         = GFM.Theta_Eg0(i_GFM);

    switch Y_DG(GFM_indices(i_GFM),3)
        case 'Droop'
        
            Control.VIM(i_GFM)          = 0; %no VIM used as Grid Angle Measurement
            Control.PLL(i_GFM)          = 0; %no PLL used as Grid Angle Measurement
            Control.MMU_QVDroop(i_GFM)  = 1; %MMU: Q-V Droop enabled 
            Control.MMU_VSM(i_GFM)      = 0; %MMU: VSM AVR not enabled
            Control.dVOC(i_GFM)         = 0; %MMU: dVOC not enabled 
            Control.SU_PfDroop(i_GFM)   = 1; %SU: P-f Droop enabled 
            Control.SU_Matching(i_GFM)  = 0; %SU: Matching not enabled 
            Control.SU_VSM(i_GFM)       = 0; %SU: VSM not enabled 
            Control.MMU(i_GFM)          = 1; %MMU: Q-V Droop selected for Eg_set
            Control.SU_theta(i_GFM)     = 2; %SU: P-f Droop's theta selected 
            Control.SU_omega(i_GFM)     = 2; %SU: P-f Droop's omega selected
            Control.omega_set(i_GFM)    = 2; %1 pu selected as the w_set
            
        case 'Droop+PLL'
        
            Control.VIM(i_GFM)          = 0; %no VIM used as Grid Angle Measurement
            Control.PLL(i_GFM)          = 1; %PLL used as Grid Angle Measurement
            Control.MMU_QVDroop(i_GFM)  = 1; %MMU: Q-V Droop enabled 
            Control.MMU_VSM(i_GFM)      = 0; %MMU: VSM AVR not enabled
            Control.dVOC(i_GFM)         = 0; %MMU: dVOC not enabled 
            Control.SU_PfDroop(i_GFM)   = 1; %SU: P-f Droop enabled 
            Control.SU_Matching(i_GFM)  = 0; %SU: Matching not enabled 
            Control.SU_VSM(i_GFM)       = 0; %SU: VSM not enabled 
            Control.MMU(i_GFM)          = 1; %MMU: Q-V Droop selected for Eg_set
            Control.SU_theta(i_GFM)     = 2; %SU: P-f Droop's theta selected 
            Control.SU_omega(i_GFM)     = 2; %SU: P-f Droop's omega selected
            Control.omega_set(i_GFM)    = 3; %w_PLL selected as the w_set
            
        case 'Droop+VIM'
        
            Control.VIM(i_GFM)          = 1; %VIM used as Grid Angle Measurement
            Control.PLL(i_GFM)          = 0; %no PLL used as Grid Angle Measurement
            Control.MMU_QVDroop(i_GFM)  = 1; %MMU: Q-V Droop enabled 
            Control.MMU_VSM(i_GFM)      = 0; %MMU: VSM AVR not enabled
            Control.dVOC(i_GFM)         = 0; %MMU: dVOC not enabled 
            Control.SU_PfDroop(i_GFM)   = 1; %SU: P-f Droop enabled 
            Control.SU_Matching(i_GFM)  = 0; %SU: Matching not enabled 
            Control.SU_VSM(i_GFM)       = 0; %SU: VSM not enabled 
            Control.MMU(i_GFM)          = 1; %MMU: Q-V Droop selected for Eg_set
            Control.SU_theta(i_GFM)     = 2; %SU: P-f Droop's theta selected 
            Control.SU_omega(i_GFM)     = 2; %SU: P-f Droop's omega selected
            Control.omega_set(i_GFM)    = 1; %w_VIM selected as the w_set
            
        case 'dVOC'
            
            Control.VIM(i_GFM)          = 0; %no VIM used as Grid Angle Measurement 
            Control.PLL(i_GFM)          = 0; %no PLL used as Grid Angle Measurement
            Control.MMU_QVDroop(i_GFM)  = 0; %MMU: Q-V Droop not enabled
            Control.MMU_VSM(i_GFM)      = 0; %MMU: VSM AVR not enabled
            Control.dVOC(i_GFM)         = 1; %MMU: dVOC enabled 
            Control.SU_PfDroop(i_GFM)   = 0; %SU: P-f Droop not enabled 
            Control.SU_Matching(i_GFM)  = 0; %SU: Matching not enabled 
            Control.SU_VSM(i_GFM)       = 0; %SU: VSM not enabled 
            Control.MMU(i_GFM)          = 2; %MMU: dVOC selected for Eg_set
            Control.SU_theta(i_GFM)     = 1; %SU: dVOC's theta selected 
            Control.SU_omega(i_GFM)     = 1; %SU: dVOC's omega selected
            Control.omega_set(i_GFM)    = 2; %1 pu selected as the w_set
            
        case 'Matching'
             
            Control.VIM(i_GFM)          = 0; %no VIM used as Grid Angle Measurement 
            Control.PLL(i_GFM)          = 0; %no PLL used as Grid Angle Measurement
            Control.MMU_QVDroop(i_GFM)  = 0; %MMU: Q-V Droop not enabled
            Control.MMU_VSM(i_GFM)      = 0; %MMU: VSM AVR not enabled
            Control.dVOC(i_GFM)         = 0; %MMU: dVOC not enabled 
            Control.SU_PfDroop(i_GFM)   = 0; %SU: P-f Droop not enabled 
            Control.SU_Matching(i_GFM)  = 1; %SU: Matching enabled 
            Control.SU_VSM(i_GFM)       = 0; %SU: VSM not enabled 
            Control.MMU(i_GFM)          = 3; %MMU: Eg_set_0 selected 
            Control.SU_theta(i_GFM)     = 3; %SU: Matching's theta selected 
            Control.SU_omega(i_GFM)     = 3; %SU: Matching's omega selected
            Control.omega_set(i_GFM)    = 2; %1 pu selected as the w_set
            
        case 'VSM'
             
            Control.VIM(i_GFM)          = 0; %no VIM used as Grid Angle Measurement
            Control.PLL(i_GFM)          = 0; %no PLL used as Grid Angle Measurement
            Control.MMU_QVDroop(i_GFM)  = 0; %MMU: Q-V Droop not enabled
            Control.MMU_VSM(i_GFM)      = 1; %MMU: VSM AVR enabled
            Control.dVOC(i_GFM)         = 0; %MMU: dVOC not enabled 
            Control.SU_PfDroop(i_GFM)   = 0; %SU: P-f Droop not enabled 
            Control.SU_Matching(i_GFM)  = 0; %SU: Matching not enabled 
            Control.SU_VSM(i_GFM)       = 1; %SU: VSM enabled 
            Control.MMU(i_GFM)          = 4; %MMU: VSM'S AVR voltage reference is selected 
            Control.SU_theta(i_GFM)     = 4; %SU: VSM's theta selected 
            Control.SU_omega(i_GFM)     = 4; %SU: VSM's omega selected
            Control.omega_set(i_GFM)    = 2; %1 pu selected as the w_set  
        
        otherwise
            Control.VIM(i_GFM)          = 0; %VIM used as Grid Angle Measurement 
            Control.PLL(i_GFM)          = 0; %No PLL used as Grid Angle Measurement 
            Control.SU_theta(i_GFM)     = 0; %SU: VIM's theta selected 
            Control.SU_omega(i_GFM)     = 0; %SU: VIM's omega selected
            
            Control.MMU_QVDroop(i_GFM)  = 0; %MMU: Q-V Droop enabled 
            Control.MMU_VSM(i_GFM)      = 0; %MMU: VSM AVR not enabled
            Control.dVOC(i_GFM)         = 0; %MMU: dVOC not enabled 
            Control.SU_PfDroop(i_GFM)   = 0; %SU: P-f Droop enabled 
            Control.SU_Matching(i_GFM)  = 0; %SU: Matching not enabled 
            Control.SU_VSM(i_GFM)       = 0; %SU: VSM not enabled 
            Control.MMU(i_GFM)          = 0; %MMU: Q-V Droop selected for Eg_set
            Control.omega_set(i_GFM)    = 0; %1 pu selected as the w_set
            
    end
end

% VSC-GFL Initializations:
GFL_indices = (find(Y_network(:,2)==3)).';
for i_GFL = 1:n_GFL
    GFL.Pref(i_GFL) = bus_sol(GFL_indices(i_GFL),4); %in pu
    GFL.Qref(i_GFL) = bus_sol(GFL_indices(i_GFL),5); %in pu
    GFL.V(i_GFL) = bus_sol(GFL_indices(i_GFL),2); %in pu
    GFL.Theta_Eg0(i_GFL) = bus_sol(GFL_indices(i_GFL),3)*pi/180;
    
    GFL.Eg0(i_GFL) = GFL.V(i_GFL)*exp(1i*GFL.Theta_Eg0(i_GFL));
    GFL.Eg0dq(i_GFL) = GFL.Eg0(i_GFL)*exp(-1i*GFL.Theta_Eg0(i_GFL));
    GFL.Ig0(i_GFL) = (GFL.Pref(i_GFL) - 1i*GFL.Qref(i_GFL))/conj(GFL.Eg0(i_GFL));
    GFL.Ig0dq(i_GFL) = GFL.Ig0(i_GFL)*exp(-1i*GFL.Theta_Eg0(i_GFL));
    GFL.Ic0(i_GFL) = GFL.Eg0(i_GFL)*(1i*GFL.C_pu);
    GFL.Ic0dq(i_GFL) = GFL.Ic0(i_GFL)*exp(-1i*GFL.Theta_Eg0(i_GFL));
    GFL.Is0(i_GFL) = GFL.Ig0(i_GFL) + GFL.Ic0(i_GFL);
    GFL.Is0dq(i_GFL) = GFL.Ig0dq(i_GFL) + GFL.Ic0dq(i_GFL);
    GFL.Vm0(i_GFL) = GFL.Eg0(i_GFL) + GFL.Is0(i_GFL)*(GFL.R1_pu+1i*GFL.L1_pu);
    GFL.Vm0dq(i_GFL) = GFL.Vm0(i_GFL)*exp(-1i*GFL.Theta_Eg0(i_GFL));
    GFL.Theta_Vm0(i_GFL) = angle(GFL.Vm0(i_GFL));
%     GFL.Ig0(i_GFL) = (GFL.Pref(i_GFL) - i*GFL.Qref(i_GFL))/conj(GFL.eg0(i_GFL));
%     GFL.It(i_GFL) = abs(GFL.Ig0(i_GFL));
%     GFL.Ic(i_GFL) = GFL.eg0(i_GFL)*(i*GFL.C_pu);
%     GFL.Phi(i_GFL) = angle(GFL.Pref(i_GFL) - i*GFL.Qref(i_GFL));
%     GFL.isdq0(i_GFL) = GFL.Ig0(i_GFL)+GFL.Ic(i_GFL);
%     GFL.vmdq0(i_GFL) = (GFL.eg0(i_GFL)+GFL.isdq0(i_GFL)*((GFL.R1_pu+i*GFL.L1_pu)));
%     GFL.delta0(i_GFL) = angle(GFL.vmdq0(i_GFL));
%     
%     GFL.vmd_0(i_GFL) = real(GFL.vmdq0(i_GFL));
%     GFL.vmq_0(i_GFL) = imag(GFL.vmdq0(i_GFL));
% 
%     GFL.isd_0(i_GFL) = real(GFL.isdq0(i_GFL));
%     GFL.isq_0(i_SG) = imag(GFL.isdq0(i_GFL));

    % Grid voltage
    GFL.Eg0_d(i_GFL)  = real(GFL.Eg0dq(i_GFL)); 
    GFL.Eg0_q(i_GFL)  = imag(GFL.Eg0dq(i_GFL));
    GFL.Eg0_dq(i_GFL,:) = [GFL.Eg0_d(i_GFL); GFL.Eg0_q(i_GFL); 0]'; %column_1 is Ed, column_2 is Eq and column_3 is E0

    GFL.P0(i_GFL) = GFL.Pref(i_GFL);
    GFL.Q0(i_GFL) = GFL.Qref(i_GFL);
    GFL.Ws0       = Parameters.wb;

    % Ig
    GFL.Ig0_d(i_GFL)  = real(GFL.Ig0dq(i_GFL));
    GFL.Ig0_q(i_GFL)  = imag(GFL.Ig0dq(i_GFL)); 
    GFL.Ig0_dq(i_GFL,:) = [GFL.Ig0_d(i_GFL); GFL.Ig0_q(i_GFL) ; 0]'; %column_1 is Igd, column_2 is Igq and column_3 is Ig0

    GFL.Is0_d(i_GFL)  = real(GFL.Is0dq(i_GFL)); 
    GFL.Is0_q(i_GFL)  = imag(GFL.Is0dq(i_GFL)); 
    GFL.Is0_dq(i_GFL,:) = [GFL.Is0_d(i_GFL); GFL.Is0_q(i_GFL); 0]';

    % Vconv
    GFL.Vm0_d(i_GFL)  = real(GFL.Vm0dq(i_GFL));
    GFL.Vm0_q(i_GFL)  = imag(GFL.Vm0dq(i_GFL));
    GFL.Vm0_dq(i_GFL,:)  = [GFL.Vm0_d(i_GFL); GFL.Vm0_q(i_GFL); 0]';
 
    GFL.theta0(i_GFL) = GFL.Theta_Eg0(i_GFL);
    GFL.M_transform_inv = [   sin(GFL.theta0(i_GFL))     ,    cos(GFL.theta0(i_GFL))     ,    1     ;
                          sin(GFL.theta0(i_GFL)-2*pi/3)     ,  cos(GFL.theta0(i_GFL)-2*pi/3)  , 1  ;
                          sin(GFL.theta0(i_GFL)+2*pi/3)     ,  cos(GFL.theta0(i_GFL)+2*pi/3)  ,  1 ];
%
    GFL.Eg0_abc(i_GFL,:) = (GFL.M_transform_inv * GFL.Eg0_dq(i_GFL,:)')';
    GFL.Vm0_abc(i_GFL,:) = (GFL.M_transform_inv * GFL.Vm0_dq(i_GFL,:)')';
    GFL.Ig0_abc(i_GFL,:) = (GFL.M_transform_inv * GFL.Ig0_dq(i_GFL,:)')';
    GFL.Is0_abc(i_GFL,:) = (GFL.M_transform_inv * GFL.Is0_dq(i_GFL,:)')';

    % Inner Current Loop Integrator Initilization: 
    GFL_Inner_current_loop.Int_ini_d(i_GFL) = GFL.Vm0_d(i_GFL) - (GFL.Eg0_d(i_GFL) .* GFL_Inner_current_loop.Kffv(i_GFL)) + (GFL.Is0_q(i_GFL) .* GFL.L1_pu * 1) ;
    GFL_Inner_current_loop.Int_ini_q(i_GFL) = GFL.Vm0_q(i_GFL) - (GFL.Eg0_q(i_GFL) .* GFL_Inner_current_loop.Kffv(i_GFL)) - (GFL.Is0_d(i_GFL) .* GFL.L1_pu * 1) ;


    % Filter Currents and voltage initialization:
    GFL_Filter.I_ini(i_GFL,:)             = GFL.Is0_abc(i_GFL,:).*GFL.Ib ;  %
    GFL_Filter.V_ini(i_GFL,:)             = GFL.Eg0_abc(i_GFL,:).*GFL.Ub ;  %
    GFL_Line.I_ini(i_GFL,:)               = GFL.Ig0_abc(i_GFL,:).*GFL.Ib ;  %
 
    % GFL initializations:
    M_pll_0(i_GFL)                  = (-1 - GFL.Eg0_q(i_GFL).*PLL.Kp)./PLL.Ki;
    M_d_0(i_GFL)                    = GFL.Is0_d(i_GFL) ./ GFL.vdc_Ki;
    M_q_0(i_GFL)                    = GFL.Is0_q(i_GFL) ./ GFL.KQ;
    
    switch Y_DG(GFL_indices(i_GFL),3)
        case 'PLL'
        
            GFL_Control.VIM(i_GFL)          = 0; %no VIM used as Grid Angle Measurement 
            GFL_Control.PLL(i_GFL)          = 1; %PLL used as Grid Angle Measurement 
            GFL_Control.SU_theta(i_GFL)     = 1; %SU: PLL's theta selected 
            GFL_Control.SU_omega(i_GFL)     = 1; %SU: PLL's omega selected
                        
        case 'VIM'
            
            GFL_Control.VIM(i_GFL)          = 1; %VIM used as Grid Angle Measurement 
            GFL_Control.PLL(i_GFL)          = 0; %No PLL used as Grid Angle Measurement 
            GFL_Control.SU_theta(i_GFL)     = 2; %SU: VIM's theta selected 
            GFL_Control.SU_omega(i_GFL)     = 2; %SU: VIM's omega selected
                        
    end
end

%% ========================================================================
%  Section VII: Initializing the Masks
%  ========================================================================
if Dynamic_sim ==1
open_system(model.name)

[nodes_Un,nodes_Type,GFM_nType,SG_nType] = nodes_NP (Y_network1,Y_DER,Y_TR,Lines_BP,IB_BP,SM_BP,VSC_BP,n_GFM,n_SG);
% The Upstream Network (Infinite Bus) Mask:
if sflag_IB == 1
    set_param([model.name,'/Upstream Network'],...
        'Voltage',  num2str(IB_BP(1)),...
        'Frequency', num2str(Parameters.fb),...
        'BaseVoltage', num2str(IB_BP(1)),...
        'BusType', 'swing');   % always cnsidered as a swing   
end

% The SGs' Masks
if n_SG > 0 
    for i_SG = 1:n_SG
        set_param([model.name,'/SM_',num2str(i_SG)],...
            'index_SM', num2str(i_SG),...
            'Parameters_fb', num2str(Parameters.fb),...
            'SM_Sn', num2str(SM.Sn),...
            'SM_Un', num2str(SM.Un),...
            'SM_wr', num2str(SM.wr),...
            'SM_Ra', num2str(SM.Ra),...
            'SM_Ll', num2str(SM.Ll),...
            'SM_Lad', num2str(SM.Lad),...
            'SM_Laq', num2str(SM.Laq),...
            'SM_L1d', num2str(SM.L1d),...
            'SM_R1d', num2str(SM.R1d),...
            'SM_L1q', num2str(SM.L1q),...
            'SM_R1q', num2str(SM.R1q),...
            'SM_L2q', num2str(SM.L2q),...
            'SM_R2q', num2str(SM.R2q),...
            'SM_Rfd', num2str(SM.Rfd),...
            'SM_Lfd', num2str(SM.Lfd),...
            'SM_H', num2str(SM.H(i_SG)),...
            'SM_K_D', num2str(SM.K_D(i_SG)),...
            'SM_Es_base', num2str(SM.Es_base),...
            'SM_Is_base', num2str(SM.Is_base),...
            'SM_Pt', num2str(SM.Pt(i_SG)),...
            'SM_Qt', num2str(SM.Qt(i_SG)),...
            'SM_Et', num2str(SM.Et(i_SG)),...
            'SM_It', num2str(SM.It(i_SG)),...
            'SM_Te_0', num2str(SM.Te_0(i_SG)),...
            'SM_alpha', num2str(SM.alpha(i_SG)),...
            'SM_Phi', num2str(SM.Phi(i_SG)),...
            'SM_delta0', num2str(SM.delta0(i_SG)),...
            'SM_theta_0', num2str(SM.theta_0(i_SG)),...
            'SM_ed_0', num2str(SM.ed_0(i_SG)),...
            'SM_eq_0', num2str(SM.eq_0(i_SG)),...
            'SM_id_0', num2str(SM.id_0(i_SG)),...
            'SM_iq_0', num2str(SM.iq_0(i_SG)),...
            'SM_efd_0', num2str(SM.efd_0(i_SG)),...
            'SM_ifd_0', num2str(SM.ifd_0(i_SG)),...
            'SM_Psi_d0', num2str(SM.Psi_d0(i_SG)),...
            'SM_Psi_q0', num2str(SM.Psi_q0(i_SG)),...
            'SM_Psi_fd0', num2str(SM.Psi_fd0(i_SG)),...
            'SM_Psi_1d0', num2str(SM.Psi_1d0(i_SG)),...
            'SM_Psi_1q0', num2str(SM.Psi_1q0(i_SG)),...
            'SM_Psi_2q0', num2str(SM.Psi_2q0(i_SG)),...
            'SM_mp', num2str(SM.mp(i_SG)),...
            'SM_Tg', num2str(SM.Tg),...
            'SM_Tr', num2str(SM.Tr),...
            'SM_Ka', num2str(SM.Ka),...
            'SM_Ta', num2str(SM.Ta),...
            'SM_Ke', num2str(SM.Ke),...
            'SM_Te', num2str(SM.Te),...
            'SM_Tb', num2str(SM.Tb),...
            'SM_Tc', num2str(SM.Tc),...
            'SM_Kf', num2str(SM.Kf),...
            'SM_Tf', num2str(SM.Tf),...
            'SM_Efmin', num2str(SM.Efmin),...
            'SM_Efmax', num2str(SM.Efmax),...
            'SM_Kp', num2str(SM.Kp),...
            'SM_T_LP', num2str(SM.T_LP),...
            'SM_K_PSS', num2str(SM.K_PSS),...
            'SM_T_HP', num2str(SM.T_HP),...
            'SM_T1n', num2str(SM.T1n),...
            'SM_T1d', num2str(SM.T1d),...
            'SM_T2n', num2str(SM.T2n),...
            'SM_T2d', num2str(SM.T2d),...
            'SM_VSmax', num2str(SM.VSmax),...
            'SM_VSmin', num2str(SM.VSmin),...
            'delta_Pref', num2str(delta_Pref(i_SG)),...
            'T_Pref', num2str(T_Pref(i_SG)),...
            'delta_Vref', num2str(delta_Vref(i_SG)),...
            'T_Vref', num2str(T_Vref(i_SG))); 
        set_param([model.name,'/SM_',num2str(i_SG),'/Synchronous Machine pu Fundamental'],...
            'BusType', SG_nType);        
    end
end

% The GFMs' Masks:
if n_GFM > 0 
    for i_GFM = 1:n_GFM
        set_param([model.name,'/GFM_',num2str(i_GFM)],...
                'index_c', num2str(i_GFM),...
                'Parameters_fb', num2str(Parameters.fb),...
                'Parameters_wb', num2str(Parameters.wb),...
                'GFM_Un', num2str(GFM.Un),...
                'GFM_Un_dc', num2str(GFM.Un_dc),...
                'GFM_Sn', num2str(GFM.Sn),...
                'GFM_Ib', num2str(GFM.Ib),...
                'GFM_Ub', num2str(GFM.Ub),...
                'GFM_Ub_dc', num2str(GFM.Ub_dc),...
                'GFM_Ib_dc', num2str(GFM.Ib_dc),...
                'GFM_Pref', num2str(GFM.Pref(i_GFM)),...
                'GFM_Qref', num2str(GFM.Qref(i_GFM)),...
                'GFM_Vref', num2str(GFM.Vref(i_GFM)),...
                'GFM_Theta_Eg0', num2str(GFM.Theta_Eg0(i_GFM)),...
                'GFM_Vdcref', num2str(GFM.Vdcref),...
                'Delta_Pref_GFM', num2str(Delta_Pref_GFM(i_GFM)),...
                'Tsw_Pref_GFM', num2str(Tsw_Pref_GFM(i_GFM)),...
                'Delta_Vref_GFM', num2str(Delta_Vref_GFM(i_GFM)),...
                'Tsw_Vref_GFM', num2str(Tsw_Vref_GFM(i_GFM)),...
                'DC_link_R', num2str(DC_link.R),...
                'DC_link_C', num2str(DC_link.C),...
                'DC_link_Vdc_init', num2str(DC_link.Vdc_init),...
                'GFM_R1_pu', num2str(GFM.R1_pu),...
                'GFM_L1_pu', num2str(GFM.L1_pu),...
                'GFM_R2_pu', num2str(GFM.R2_pu),...
                'GFM_L2_pu', num2str(GFM.L2_pu),...
                'GFM_C_pu', num2str(GFM.C_pu),...
                'GFM_R1', num2str(GFM.R1),...
                'GFM_L1', num2str(GFM.L1),...
                'GFM_R2', num2str(GFM.R2),...
                'GFM_L2', num2str(GFM.L2),...
                'GFM_C', num2str(GFM.C),...
                'Controllers', num2str(Y_DG(GFM_indices(i_GFM),3)),....
                'Control_VIM', num2str(Control.VIM(i_GFM)),...
                'Control_PLL', num2str(Control.PLL(i_GFM)),...
                'Control_MMU_QVDroop', num2str(Control.MMU_QVDroop(i_GFM)),...
                'Control_MMU_VSM', num2str(Control.MMU_VSM(i_GFM)),...
                'Control_dVOC', num2str(Control.dVOC(i_GFM)),...
                'Control_MMU', num2str(Control.MMU(i_GFM)),...
                'Control_SU_PfDroop', num2str(Control.SU_PfDroop(i_GFM)),...
                'Control_SU_Matching', num2str(Control.SU_Matching(i_GFM)),...
                'Control_SU_VSM', num2str(Control.SU_VSM(i_GFM)),...
                'Control_SU_theta', num2str(Control.SU_theta(i_GFM)),...
                'Control_SU_omega', num2str(Control.SU_omega(i_GFM)),...
                'Control_omega_set', num2str(Control.omega_set(i_GFM)),...
                'Droop_mp', num2str(Droop.mp(i_GFM)),...
                'Droop_nq', num2str(Droop.nq(i_GFM)),...
                'Droop_wc', num2str(Droop.wc(i_GFM)),...
                'Droop_Int_ini', num2str(Droop.Int_ini),...
                'q_filter_wf', num2str(q_filter.wf(i_GFM)),...
                'q_filter_Int_ini', num2str(q_filter.Int_ini(i_GFM)),...
                'p_filter_Int_ini', num2str(p_filter.Int_ini(i_GFM)),...
                'Vdc_loop_Ki', num2str(Vdc_loop.Ki),...
                'Vdc_loop_Kp', num2str(Vdc_loop.Kp),...
                'Vdc_loop_Kffi', num2str(Vdc_loop.Kffi(i_GFM)),...
                'Vdc_loop_Int_ini', num2str(Vdc_loop.Int_ini),...
                'Damping_resistor_Enable', num2str(Damping_resistor.Enable(i_GFM)),...
                'Damping_resistor_Rv', num2str(Damping_resistor.Rv),...
                'Damping_resistor_W_RV', num2str(Damping_resistor.W_RV(i_GFM)),...
                'VI_Enable', num2str(VI.Enable(i_GFM)),...
                'VI_Rv_pu', num2str(VI.Rv_pu(i_GFM)),...
                'VI_Lv_pu', num2str(VI.Lv_pu(i_GFM)),...
                'TVI_K_RV', num2str(TVI.K_RV(i_GFM)),...
                'TVI_Rho_XR', num2str(TVI.Rho_XR(i_GFM)),...
                'TVI_I_thresh_pu', num2str(TVI.I_thresh_pu(i_GFM)),...
                'Inner_voltage_loop_Kp', num2str(Inner_voltage_loop.Kp),...
                'Inner_voltage_loop_Ki', num2str(Inner_voltage_loop.Ki),...
                'Inner_voltage_loop_Kffi', num2str(Inner_voltage_loop.Kffi(i_GFM)),...
                'Inner_current_loop_Kp', num2str(Inner_current_loop.Kp),...
                'Inner_current_loop_Ki', num2str(Inner_current_loop.Ki),...
                'Inner_current_loop_Kffv', num2str(Inner_current_loop.Kffv(i_GFM)),...
                'SFR_Enable', num2str(SFR.Enable(i_GFM)),...
                'SFR_Rc', num2str(SFR.Rc(i_GFM)),...
                'SFR_dP_max', num2str(SFR.dP_max(i_GFM)),...
                'GFM_theta', num2str(theta0(i_GFM)),...
                'Inner_voltage_loop_Int_ini_d', num2str(Inner_voltage_loop.Int_ini_d(i_GFM)),...
                'Inner_voltage_loop_Int_ini_q', num2str(Inner_voltage_loop.Int_ini_q(i_GFM)),...
                'Inner_current_loop_Int_ini_d', num2str(Inner_current_loop.Int_ini_d(i_GFM)),...
                'Inner_current_loop_Int_ini_q', num2str(Inner_current_loop.Int_ini_q(i_GFM)),...
                'Damping_resistor_Int_ini', num2str(Damping_resistor.Int_ini(i_GFM)),...
                'Filter_I_ini', ['[ ',num2str(Filter.I_ini(i_GFM,:)),' ]'],...
                'Filter_V_ini', ['[ ',num2str(Filter.V_ini(i_GFM,:)),' ]'],...
                'Line_I_ini', ['[ ',num2str(Line.I_ini(i_GFM,:)),' ]'],...
                'Is0_dq0', ['[ ',num2str(Is0_dq(i_GFM,:)),' ]'],...
                'Eg0_dq0', ['[ ',num2str(Eg0_dq(i_GFM,:)),' ]'],...
                'Ig0_dq0', ['[ ',num2str(Ig0_dq(i_GFM,:)),' ]'],...
                'Vm0_dq0', ['[ ',num2str(Vm0_dq(i_GFM,:)),' ]'],...
                'dVOC_eta', num2str(dVOC.eta(i_GFM)),...
                'dVOC_alpha', num2str(dVOC.alpha(i_GFM)),...
                'Matching_K_theta', num2str(Matching.K_theta(i_GFM)),...               
                'VSM_D', num2str(VSM.D(i_GFM)),...
                'VSM_J', num2str(VSM.J(i_GFM)),...
                'VSM_Dq', num2str(VSM.Dq(i_GFM)),...
                'VSM_K', num2str(VSM.K(i_GFM)),...
                'dVOC_ini_theta', num2str(dVOC.ini_theta(i_GFM)),...
                'Matching_Int_ini', num2str(Matching.Int_ini(i_GFM)));
            set_param([model.name,'/Init_GFM_',num2str(i_GFM)],...
                'Voltage',  num2str(GFM.Un),...
                'Frequency', num2str(Parameters.fb),...
                'BaseVoltage', num2str(GFM.Un),...
                'BusType', GFM_nType,...                
                'Pref', num2str(GFM.Pref(i_GFM)*GFM.Sn),...
                'Qref', num2str(GFM.Qref(i_GFM)*GFM.Sn));             
    end
end

% The GFLs' Masks:
if n_GFL > 0 
    for i_GFL = 1:n_GFL
        set_param([model.name,'/GFL_',num2str(i_GFL)],...
                'index_c', num2str(i_GFL),...
                'Parameters_fb', num2str(Parameters.fb),...
                'Parameters_wb', num2str(Parameters.wb),...
                'GFL_Un', num2str(GFL.Un),...
                'GFL_Un_dc', num2str(GFL.Un_dc),...
                'GFL_Sn', num2str(GFL.Sn),...
                'GFL_Ib', num2str(GFL.Ib),...
                'GFL_Ub', num2str(GFL.Ub),...
                'GFL_Ub_dc', num2str(GFL.Ub_dc),...
                'GFL_Ib_dc', num2str(GFL.Ib_dc),...
                'GFL_Pref', num2str(GFL.Pref(i_GFL)),...
                'GFL_Qref', num2str(GFL.Qref(i_GFL)),...
                'GFL_V', num2str(GFL.V(i_GFL)),...
                'GFL_Theta_Eg0', num2str(GFL.Theta_Eg0(i_GFL)),...
                'GFL_Vdcref', num2str(GFL.Vdcref),...
                'Delta_Pref_GFL', num2str(Delta_Pref_GFL(i_GFL)),...
                'Tsw_Pref_GFL', num2str(Tsw_Pref_GFL(i_GFL)),...
                'Delta_Qref_GFL', num2str(Delta_Qref_GFL(i_GFL)),...
                'Tsw_Qref_GFL', num2str(Tsw_Qref_GFL(i_GFL)),...
                'GFL_DC_link_R', num2str(GFL_DC_link.R),...
                'GFL_DC_link_C', num2str(GFL_DC_link.C),...
                'GFL_DC_link_Vdc_init', num2str(GFL_DC_link.Vdc_init),...
                'GFL_R1_pu', num2str(GFL.R1_pu),...
                'GFL_L1_pu', num2str(GFL.L1_pu),...
                'GFL_R2_pu', num2str(GFL.R2_pu),...
                'GFL_L2_pu', num2str(GFL.L2_pu),...
                'GFL_C_pu', num2str(GFL.C_pu),...
                'GFL_R1', num2str(GFL.R1),...
                'GFL_L1', num2str(GFL.L1),...
                'GFL_R2', num2str(GFL.R2),...
                'GFL_L2', num2str(GFL.L2),...
                'GFL_C', num2str(GFL.C),...
                'Controllers', num2str(Y_DG(GFL_indices(i_GFL),3)),....
                'GFL_Control_VIM', num2str(GFL_Control.VIM(i_GFL)),...
                'GFL_Control_PLL', num2str(GFL_Control.PLL(i_GFL)),...
                'GFL_Control_SU_theta', num2str(GFL_Control.SU_theta(i_GFL)),...
                'GFL_Control_SU_omega', num2str(GFL_Control.SU_omega(i_GFL)),...
                'GFL_vdc_Kp', num2str(GFL.vdc_Kp),...
                'GFL_vdc_Ki', num2str(GFL.vdc_Ki),...
                'GFL_KQ', num2str(GFL.KQ),...
                'PLL_Kp', num2str(PLL.Kp),...
                'PLL_Ki', num2str(PLL.Ki),...
                'M_pll_0', num2str(M_pll_0),...
                'VIM_wc', num2str(VIM.wc),...
                'VIM_wc_ini', num2str(VIM.wc_ini),...
                'VIM_Iqd_sat', num2str(VIM.Iqd_sat),...
                'VIM_Rr_pu', num2str(VIM.Rr_pu),...
                'VIM_Lr_pu', num2str(VIM.Lr_pu),...
                'VIM_Lm_pu', num2str(VIM.Lm_pu),...
                'VIM_Rr', num2str(VIM.Rr),...
                'VIM_Lr', num2str(VIM.Lr),...
                'VIM_Lm', num2str(VIM.Lm),...
                'VIM_Kp', num2str(VIM.Kp),...
                'VIM_Kd', num2str(VIM.Kd),...
                'VIM_Nd', num2str(VIM.Nd),...
                'VIM_H', num2str(VIM.H),...
                'VIM_J', num2str(VIM.J),...
                'VIM_D', num2str(VIM.D),...
                'VIM_Inertia_ini', num2str(VIM.Inertia_ini),...
                'VIM_Ke_ini', num2str(VIM.Ke_ini),...
                'VIM_T_sync', num2str(VIM.T_sync),...
                'VIM_f0', num2str(VIM.f0),...
                'VIM_w0', num2str(VIM.w0),... 
                'GFL_Inner_current_loop_Kp', num2str(GFL_Inner_current_loop.Kp),...
                'GFL_Inner_current_loop_Ki', num2str(GFL_Inner_current_loop.Ki),...
                'GFL_Inner_current_loop_Kffv', num2str(GFL_Inner_current_loop.Kffv(i_GFL)),...
                'GFL_theta', num2str(GFL.theta0(i_GFL)),...
                'GFL_Inner_current_loop_Int_ini_d', num2str(GFL_Inner_current_loop.Int_ini_d(i_GFL)),...
                'GFL_Inner_current_loop_Int_ini_q', num2str(GFL_Inner_current_loop.Int_ini_q(i_GFL)),...
                'GFL_Filter_I_ini', ['[ ',num2str(GFL_Filter.I_ini(i_GFL,:)),' ]'],...
                'GFL_Filter_V_ini', ['[ ',num2str(GFL_Filter.V_ini(i_GFL,:)),' ]'],...
                'GFL_Line_I_ini', ['[ ',num2str(GFL_Line.I_ini(i_GFL,:)),' ]'],...
                'GFL_Is0_dq0', ['[ ',num2str(GFL.Is0_dq(i_GFL,:)),' ]'],...
                'GFL_Eg0_dq0', ['[ ',num2str(GFL.Eg0_dq(i_GFL,:)),' ]'],...
                'GFL_Ig0_dq0', ['[ ',num2str(GFL.Ig0_dq(i_GFL,:)),' ]'],...
                'GFL_Vm0_dq0', ['[ ',num2str(GFL.Vm0_dq(i_GFL,:)),' ]'],...
                'M_d_0', num2str(M_d_0(i_GFL)),...
                'M_q_0', num2str(M_q_0(i_GFL)));
            set_param([model.name,'/Init_GFL_',num2str(i_GFL)],...
                'Voltage',  num2str(GFL.Un),...
                'Frequency', num2str(Parameters.fb),...
                'BaseVoltage', num2str(GFL.Un),...
                'BusType', "PQ",...
                'Pref', num2str(GFL.Pref(i_GFL)*GFL.Sn),...
                'Qref', num2str(GFL.Qref(i_GFL)*GFL.Sn));               
    end
end

% The Lines' Masks
if n_lines > 0 
    for i_line = 1:n_lines
        set_param([model.name,'/Line_',num2str(i_line)],...
            'Length', num2str(Y_line(i_line,7)),...
            'Frequency', num2str(Parameters.fb),...
            'Resistances',  ['[ ',num2str(Line.R(i_line)),' ' , num2str(0.3864),' ]'],...
            'Inductances', ['[ ',num2str(Line.L(i_line)),' ', num2str(4.1264e-3),' ]'],...
            'Capacitances', ['[ ',num2str(Line.C(i_line)),' ', num2str(7.751e-9),' ]']);            
    end
end

% The Transformers' Masks
if n_TR > 0 
    for i_TR = 1:n_TR
        set_param([model.name,'/Transformer_',num2str(i_TR)],...
            'NominalPower',  ['[ ',num2str(TR.Sn(i_TR)),' ' , num2str(Parameters.fb),' ]'],...
            'Winding1', ['[ ',num2str(TR.Un1),' ', num2str(TR.R1(i_TR)),' ', num2str(TR.L1(i_TR)),' ]'],...
            'Winding2', ['[ ',num2str(TR.Un2),' ', num2str(TR.R2(i_TR)),' ', num2str(TR.L2(i_TR)),' ]'],...
            'Rm', num2str(TR.Rm),...
            'Lm', num2str(TR.Lm));            
    end
end

% The Loads' Masks
loadsim = 0;
if n_loads > 0 
    for i_node = 1:n_nodes
        if Y_network(i_node,3) == 1
            loadsim = loadsim + 1;
            set_param([model.name,'/Load_',num2str(loadsim),'/Series RLC Load'],...
                'NominalVoltage',  num2str(Load.Un),...
                'NominalFrequency', num2str(Parameters.fb),...
                'ActivePower', num2str(Load.P(i_node)),...
                'InductivePower', num2str(Load.Q(i_node))); 
            set_param([model.name,'/Load_',num2str(loadsim),'/Dynamic Load'],...
                'NominalVoltage',  ['[ ',num2str(Load.Un),' ' , num2str(Parameters.fb),' ]']); 
            set_param([model.name,'/Load_',num2str(loadsim),'/Step_P'],...
                'StepV', num2str(delta_Pload(loadsim)),...
                'StepT', num2str(T_Pload(loadsim))); 
            set_param([model.name,'/Load_',num2str(loadsim),'/Step_Q'],...
                'StepV', num2str(delta_Qload(loadsim)),...
                'StepT', num2str(T_Qload(loadsim)));                
        else
            continue;
        end
    end
end

% The Buses' Masks
if n_nodes > 0 
    for i_node = 1:n_nodes
        set_param([model.name,'/Bus_',num2str(i_node)],...
            'LabelV', ['V_Bus', num2str(i_node)],...
            'LabelI', ['I_Bus', num2str(i_node)]); 
        set_param([model.name,'/Load Flow Bus',num2str(i_node)],...
            'ID', ['BUS_', num2str(i_node)],...
            'Vbase', num2str(nodes_Un(i_node)),...
            'Vref', num2str(Y_network(i_node,5)),...
            'Vangle', num2str(Y_DER(1,6)));         
    end
end
%% ========================================================================
%  Section VIII: Running Load Flow in Simulink and applying it to
%  initialize the Simulink models
%  ========================================================================
%% Routine to run loadflow in simulink:
if reinit_LF == 1  % in case you chose to re-initialise the Loadflow
    % Step 1: Comment/Uncomment blocks
    for i_GFM = 1:n_GFM
        blockToComment1 = [model.name,'/GFM_',num2str(i_GFM)]; % Replace 'block_path_to_comment' with the path to the block you want to comment
        blockToComment2 = [model.name,'/GFM',num2str(i_GFM),'_meas']; % Replace 'block_path_to_comment' with the path to the block you want to comment
        blockToUncomment = [model.name,'/Init_GFM_',num2str(i_GFM)]; % Replace 'block_path_to_uncomment' with the path to the block you want to uncomment

        set_param(blockToComment1, 'Commented', 'on');
        set_param(blockToComment2, 'Commented', 'on');
        set_param(blockToUncomment, 'Commented', 'off');
    end
    for i_GFL = 1:n_GFL
        blockToComment1 = [model.name,'/GFL_',num2str(i_GFL)]; % Replace 'block_path_to_comment' with the path to the block you want to comment
        blockToComment2 = [model.name,'/GFL',num2str(i_GFL),'_meas']; % Replace 'block_path_to_comment' with the path to the block you want to comment
        blockToUncomment = [model.name,'/Init_GFL_',num2str(i_GFL)]; % Replace 'block_path_to_uncomment' with the path to the block you want to uncomment

        set_param(blockToComment1, 'Commented', 'on');
        set_param(blockToComment2, 'Commented', 'on');    
        set_param(blockToUncomment, 'Commented', 'off');
    end
    % Step 2: Run load flow and applies it to the model
    LF = power_loadflow('-v2',model.name,'solve');

% Step 4: Undo changes in Step 1
    for i_GFM = 1:n_GFM
        blockToComment1 = [model.name,'/GFM_',num2str(i_GFM)];
        blockToComment2 = [model.name,'/GFM',num2str(i_GFM),'_meas'];
        blockToUncomment = [model.name,'/Init_GFM_',num2str(i_GFM)];

        set_param(blockToComment1, 'Commented', 'off');
        set_param(blockToComment2, 'Commented', 'off');
        set_param(blockToUncomment, 'Commented', 'on');
    end
    for i_GFL = 1:n_GFL
        blockToComment1 = [model.name,'/GFL_',num2str(i_GFL)];
        blockToComment2 = [model.name,'/GFL',num2str(i_GFL),'_meas'];
        blockToUncomment = [model.name,'/Init_GFL_',num2str(i_GFL)];

        set_param(blockToComment1, 'Commented', 'off');
        set_param(blockToComment2, 'Commented', 'off');
        set_param(blockToUncomment, 'Commented', 'on');
    end
%    pause(30); 
end
%%
if runSim == 1
    sim(model.name);
end
% save(sprintf('SIM_CIGRE_Unstable_07.mat'),'logsout');
end
%% ========================================================================
%  Section IX: Linearization Procedure
%  ========================================================================
if ss_analysis == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user specified inputs, The code below shouldn't be modified unless
% you know what you are doing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calling the individual DG state spaces:

if sflag_IB == 1   % The slack DG is an IB
    [A_s0,B_s0,C_s0,D_s0,states_s,inputs_s,outputs_s,nStates_s,nIP_s,nOP_s] = symIB();
elseif sflag_GFM == 1   % The slack DG is a GFM
    [A_s0,B_s0,C_s0,D_s0,states_s,inputs_s,outputs_s,nStates_s,nIP_s,nOP_s] =symGFM_types(GFM_P_control,1);  % in case GFM with droop is used at the slack bus. 
elseif sflag_SG == 1   % The slack DG is a SG
    [A_s0,B_s0,C_s0,D_s0,states_s,inputs_s,outputs_s,nStates_s,nIP_s,nOP_s] = symSG(SM_PSS,1);
end

[A_Line0,B_Line0,C_Line0,D_Line0,states_Line,inputs_Line,outputs_Line,nStates_Line,nIP_Line,nOP_Line] = symLine();
[A_Load0,B_Load0,C_Load0,D_Load0,states_Load,inputs_Load,outputs_Load,nStates_Load,nIP_Load,nOP_Load] = symLoad();
[A_Node0,B_Node0,C_Node0,D_Node0,states_Node,inputs_Node,outputs_Node,nStates_Node,nIP_Node,nOP_Node] = symNode();

if n_GFL ~= 0 
    [A_GFL0,B_GFL0,C_GFL0,D_GFL0,states_GFL,inputs_GFL,outputs_GFL,nStates_GFL,nIP_GFL,nOP_GFL] = symGFL();
else
    states_GFL       = strings;
    inputs_GFL       = strings;
    outputs_GFL      = strings;
    nStates_GFL      = 0;
    nIP_GFL          = zeros(3,1);
    nOP_GFL          = zeros(3,1);
end
if (n_GFM - sflag_GFM) ~= 0 
    [A_GFM0,B_GFM0,C_GFM0,D_GFM0,states_GFM,inputs_GFM,outputs_GFM,nStates_GFM,nIP_GFM,nOP_GFM] = symGFM_types(GFM_P_control,2);
else
    states_GFM        = strings;
    inputs_GFM        = strings;
    outputs_GFM       = strings;    
    nStates_GFM       = 0;
    nIP_GFM           = zeros(3,1);
    nOP_GFM           = zeros(3,1);
end
if (n_SG - sflag_SG) ~= 0 
    [A_SG0,B_SG0,C_SG0,D_SG0,states_SG,inputs_SG,outputs_SG,nStates_SG,nIP_SG,nOP_SG] = symSG(SM_PSS,2);
else
    states_SG         = strings;
    inputs_SG         = strings;
    outputs_SG        = strings;    
    nStates_SG        = 0;
    nIP_SG            = zeros(3,1);
    nOP_SG            = zeros(3,1);
end

%% DG Parameters and Data:
[Vab] = alpha_beta(bus_sol); %computing all node voltages in the stationary reference frame alpha-beta (real-imaginary)
% Computing the initial rotation angle of the global (reference) rotating dq-reference frame (theta_g_0):
if sflag_IB == 1  % infinite bus is the slack DG
    alpha_s   = bus_sol(n_nodes-n_DG+1,3)*pi/180;
    theta_g_0 = alpha_s;
elseif sflag_GFM == 1  % A GFM is the slack DG
    alpha_s   = bus_sol(n_nodes-n_DG+1,3)*pi/180;
    theta_g_0 = alpha_s;
elseif sflag_SG == 1  % A SG is the slack DG
    E_s       = bus_sol(n_nodes-n_DG+1,2);
    alpha_s   = bus_sol(n_nodes-n_DG+1,3)*pi/180;
    Es_ab     = E_s*exp(1i*alpha_s);
    Pt_s      = bus_sol(n_nodes-n_DG+1,4);
    Qt_s      = bus_sol(n_nodes-n_DG+1,5);
    Is_ab     = (Pt_s - 1i*Qt_s)/conj(Es_ab);
    delta_s0  = angle(Es_ab + Is_ab*(SM.Ra + 1i*SM.Lq));
    theta_g_0 = delta_s0 - pi/2;
end

% IB Data:
%                      Rup(pu)       Lup(pu)        theta_g_0(rad)         wg_0(pu)   Vup_0(pu)
IB_Data             = [Y_TR(1,4)     Y_TR(1,5)      theta_g_0              1          bus_sol(n_nodes-n_DG+1,2)  Parameters.wb];

% GFM Data:
%                      1)Rf(pu)   2)Lf(pu)   3)Cf(pu)   4)Rt(pu)   5)Lt(pu)  
%                      6)wb(rad/s)  7)Cdc(pu)  8)Gdc(pu) 9)Kpdc   10)Tdc(s)  
%                      11)wc(rad/s)  12)wf(rad/s)  13)mp  14)nq  15)Kffi  
%                      16)Kffv  17)KiVL  18)KpVL  
%                      19)KiCL 20)KpCL  21)p_ref_0  22)q_ref_0  
%                      23)ve_ref_0 24)theta_0 25)Eg0_d 26)Eg0_q 27)Ig0_d 28)Ig0_q
%                      29)Is0_d 30)Is0_q 31)Vm0_d 32)Vm0_q 33)theta_g_0 
%                      34)dVOC.eta 35)dVOC.alpha 36)VSM.D 37)VSM.J
%                      38)VSM.Dq 39)VSM.K 40)Matching.K_theta
GFM_Data            = [];
for i = 1:n_GFM
GFM_Data(i,:)       = [GFM.R1_pu GFM.L1_pu GFM.C_pu Y_TR(1,4) Y_TR(1,5) ... 
                       Parameters.wb DC_link.C_pu DC_link.G_pu Vdc_loop.Kp DC_link.T_dc(i) ...
                       Droop.wc(i) q_filter.wf(i) Droop.mp(i) Droop.nq(i) Inner_voltage_loop.Kffi(i)...
                       Inner_current_loop.Kffv(i) Inner_voltage_loop.Ki Inner_voltage_loop.Kp ...
                       Inner_current_loop.Ki Inner_current_loop.Kp GFM.Pref(i) GFM.Qref(i)...
                       GFM.Vref(i) GFM.Theta_Eg0(i) Eg0_d(i) Eg0_q(i) Ig0_d(i) Ig0_q(i)  ...
                       Is0_d(i) Is0_q(i) Vm0_d(i) Vm0_q(i) theta_g_0 dVOC.eta(i) dVOC.alpha(i) ...
                       VSM.D(i) VSM.J(i) VSM.Dq(i) VSM.K(i) Matching.K_theta(i)];
end
% GFL Data:
%                      1)Rf(pu)   2)Lf(pu)   3)Cf(pu)   4)Rt(pu)   5)Lt(pu)  
%                      6)wb(rad/s)  7)Cdc(pu)  8)Gdc(pu)   9)Tdc(s)  
%                      10)Kpd  11)Kid  12)Kiq  13)Kffv  14)KpPLL  15)KiPLL  
%                      16)KiCL 17)KpCL  18)p_ref_0  19)q_ref_0  
%                      20)theta_0 21)Eg0_d 22)Eg0_q 23)Ig0_d 24)Ig0_q
%                      25)Is0_d 26)Is0_q 27)Vm0_d 28)Vm0_q 29)M_pll_0
%                      30)M_d_0 31)M_q_0 32)M_CLd_0 33)M_CLq_0
GFL_Data            = [];
for i = 1:n_GFL
GFL_Data(i,:)       = [GFL.R1_pu GFL.L1_pu GFL.C_pu Y_TR(1,4) Y_TR(1,5) ... 
                       Parameters.wb GFL_DC_link.C_pu GFL_DC_link.G_pu GFL_DC_link.T_dc(i) ...
                       GFL.vdc_Kp GFL.vdc_Ki GFL.KQ GFL_Inner_current_loop.Kffv(i) PLL.Kp PLL.Ki...
                       GFL_Inner_current_loop.Ki GFL_Inner_current_loop.Kp GFL.Pref(i) GFL.Qref(i)...
                       GFL.Theta_Eg0(i) GFL.Eg0_d(i) GFL.Eg0_q(i) GFL.Ig0_d(i) GFL.Ig0_q(i)  ...
                       GFL.Is0_d(i) GFL.Is0_q(i) GFL.Vm0_d(i) GFL.Vm0_q(i) M_pll_0(i) M_d_0(i) M_q_0(i) ...
                       GFL_Inner_current_loop.Int_ini_d(i) GFL_Inner_current_loop.Int_ini_q(i) theta_g_0];
end
% SG Data:
%                      1)wb   2)wff(pu)   3)Ra(pu)   4)Ll(pu) 5)Cg(pu) 6)Rt(pu)  7)Lt(pu)  8)Lad(pu)  9)Laq(pu)  
%                      10)Lf1d(pu)  11)Rfd(pu)  12)L11d(pu)  13)L11q(pu)  14)L22q(pu)  15)R1d(pu) 
%                      16)R1q(s)   17)R2q(s)   18)H(s) 19)KD  20)mp  
%                      21)TG  22)T_LP  23)T_HP  24)K_PSS  25)T1n  26)T1d   
%                      27)T2n  28)T2d  29)Tr  30)Ka  31)Ta  32)Ke  33)Te 34)Kfd  35)Tfd 
%                      36)Pt  37)Qt  38)Et  39)theta0 40)delta0  41)ed_0  42)eq_0
%                      43)igd_0 44)igq_0 45)Psi_d0 46)Psi_q0 47)Psi_fd0
%                      48)Psi_1d0 49)Psi_1q0 50)Psi_2q0 51)Lffd
SG_Data            = [];
for i = 1:n_SG
SG_Data(i,:)       = [Parameters.wb 1 SM.Ra SM.Ll SM.Cg_pu Y_TR(1,4) Y_TR(1,5) SM.Lads(i) SM.Laqs(i)  ... 
                      SM.Lf1d SM.Rfd SM.L11d SM.L11q SM.L22q SM.R1d ...
                      SM.R1q SM.R2q SM.H(i) SM.K_D(i)  SM.mp(i) SM.Tg ...
                      SM.T_LP SM.T_HP SM.K_PSS SM.T1n SM.T1d SM.T2n SM.T2d ...
                      SM.Tr SM.Ka SM.Ta SM.Ke SM.Te SM.Kf SM.Tf...
                      SM.Pt(i) SM.Qt(i) SM.Et(i) SM.theta_0(i) SM.delta0(i) SM.ed_0(i) SM.eq_0(i)... 
                      SM.igd_0(i) SM.igq_0(i) SM.Psi_d0(i) SM.Psi_q0(i) SM.Psi_fd0(i) ...
                      SM.Psi_1d0(i) SM.Psi_1q0(i) SM.Psi_2q0(i) SM.Lffd SM.RL_pu SM.id_0(i) SM.iq_0(i) theta_g_0];
end

%% Initial states:
init_states_IB = [1 0 theta_g_0];

init_states_SG = [];
for i = 1:n_SG
    switch SM_PSS
        case 'PSS_on'
            init_states_SG(i,:) = [SM.igd_0(i) SM.igq_0(i) SM.Psi_d0(i) SM.Psi_q0(i) ...
            SM.Psi_fd0(i) SM.Psi_1d0(i) SM.Psi_1q0(i) SM.Psi_2q0(i) 0 SM.theta_0(i) ... 
            SM.Pt(i) SM.Et(i) 0 SM.efd_0(i) 0 0 0 0 0];
        case 'PSS_off'
            init_states_SG(i,:) = [SM.igd_0(i) SM.igq_0(i) SM.Psi_d0(i) SM.Psi_q0(i) ...
            SM.Psi_fd0(i) SM.Psi_1d0(i) SM.Psi_1q0(i) SM.Psi_2q0(i) 0 SM.theta_0(i) ... 
            SM.Pt(i) 0 0 0 0];
    end
end

init_states_GFM = [];
for i = 1:n_GFM
    switch GFM_P_control
        case 'Droop'
            init_states_GFM(i,:) = [Is0_d(i) Is0_q(i) Ig0_d(i) Ig0_q(i) Eg0_d(i) Eg0_q(i) ...
                                   1 GFM.Pref(i)/1 GFM.Pref(i) GFM.Theta_Eg0(i) GFM.Qref(i) ...
                                   Inner_voltage_loop.Int_ini_d(i) Inner_voltage_loop.Int_ini_q(i)...
                                   Inner_current_loop.Int_ini_d(i) Inner_current_loop.Int_ini_q(i)];
        case 'Droop+filter'
            init_states_GFM(i,:) = [Is0_d(i) Is0_q(i) Ig0_d(i) Ig0_q(i) Eg0_d(i) Eg0_q(i) ...
                                   1 GFM.Pref(i)/1 GFM.Pref(i) 0 GFM.Theta_Eg0(i) GFM.Qref(i) ...
                                   Inner_voltage_loop.Int_ini_d(i) Inner_voltage_loop.Int_ini_q(i)...
                                   Inner_current_loop.Int_ini_d(i) Inner_current_loop.Int_ini_q(i)];            
        case 'Droop+PLL'
            init_states_GFM(i,:) = [Is0_d(i) Is0_q(i) Ig0_d(i) Ig0_q(i) Eg0_d(i) Eg0_q(i) ...
                                   1 GFM.Pref(i)/1 GFM.Pref(i) 0 GFM.Theta_Eg0(i) GFM.Qref(i) ...
                                   Inner_voltage_loop.Int_ini_d(i) Inner_voltage_loop.Int_ini_q(i)...
                                   Inner_current_loop.Int_ini_d(i) Inner_current_loop.Int_ini_q(i)];
        case 'dVOC'
            init_states_GFM(i,:) = [Is0_d(i) Is0_q(i) Ig0_d(i) Ig0_q(i) Eg0_d(i) Eg0_q(i) ...
                                   1 GFM.Pref(i)/1 GFM.Pref(i) GFM.Theta_Eg0(i) GFM.Qref(i) GFM.Vref(i) ...
                                   Inner_voltage_loop.Int_ini_d(i) Inner_voltage_loop.Int_ini_q(i)...
                                   Inner_current_loop.Int_ini_d(i) Inner_current_loop.Int_ini_q(i)];           
        case 'VSM'
            init_states_GFM(i,:) = [Is0_d(i) Is0_q(i) Ig0_d(i) Ig0_q(i) Eg0_d(i) Eg0_q(i) ...
                                   1 GFM.Pref(i)/1 0 GFM.Theta_Eg0(i) 1 ...
                                   Inner_voltage_loop.Int_ini_d(i) Inner_voltage_loop.Int_ini_q(i)...
                                   Inner_current_loop.Int_ini_d(i) Inner_current_loop.Int_ini_q(i)];            
        case 'Matching'
            init_states_GFM(i,:) = [Is0_d(i) Is0_q(i) Ig0_d(i) Ig0_q(i) Eg0_d(i) Eg0_q(i) ...
                                   1 GFM.Pref(i)/1 1 GFM.Theta_Eg0(i) ...
                                   Inner_voltage_loop.Int_ini_d(i) Inner_voltage_loop.Int_ini_q(i)...
                                   Inner_current_loop.Int_ini_d(i) Inner_current_loop.Int_ini_q(i)];            
    end
end

init_states_GFL = [];
for i = 1:n_GFL
    init_states_GFL(i,:) = [GFL.Is0_d(i) GFL.Is0_q(i) GFL.Ig0_d(i) GFL.Ig0_q(i) GFL.Eg0_d(i) GFL.Eg0_q(i) ...
                           1 GFL.Pref(i)/1 M_d_0(i) M_q_0(i) ...
                           GFL_Inner_current_loop.Int_ini_d(i) GFL_Inner_current_loop.Int_ini_q(i)...
                           M_pll_0(i) GFL.Theta_Eg0(i)];
end
%% Substituting in the individual state spaces with the init. conditions:
%%%%%%% 1. The slack DG (IB,GFM or SG): %%%%%%%%
A_s_0 = [];
B_s_0 = [];
C_s_0 = [];
D_s_0 = [];
if sflag_IB == 1
    i_DG = 1;
    [AIB0,BIB0,CIB0,DIB0] = IB_subs(Y_DER,n_nodes-n_DG+1,bus_sol,IB_Data,i_DG,A_s0,B_s0,C_s0,D_s0);
    A_s_0{1} = AIB0;
    B_s_0{1} = BIB0;
    C_s_0{1} = CIB0;
    D_s_0{1} = DIB0;
elseif sflag_GFM == 1
    i_GFM = 1;
    i_DG = 1;
    [AGFM0,BGFM0,CGFM0,DGFM0] = GFM_subs(Y_DER,bus_sol,GFM_Data,i_GFM,i_DG,A_s0,B_s0,C_s0,D_s0);
    A_s_0{1} = AGFM0;
    B_s_0{1} = BGFM0;
    C_s_0{1} = CGFM0;
    D_s_0{1} = DGFM0;
elseif sflag_SG == 1
    i_SG = 1;
    i_DG = 1;
    [ASG0,BSG0,CSG0,DSG0] = SG_subs(Y_DER,bus_sol,SG_Data,i_SG,i_DG,A_s0,B_s0,C_s0,D_s0);
    A_s_0{1} = ASG0;
    B_s_0{1} = BSG0;
    C_s_0{1} = CSG0;
    D_s_0{1} = DSG0;
end

%%%%%%%% Substituting in the different DG Types %%%%%%%%%
i_GFM = sflag_GFM; i_GFL = 0; i_SG = sflag_SG;
A_DG_0 = [];
B_DG_0 = [];
C_DG_0 = [];
D_DG_0 = [];
for i_DG = 2:n_DG
    if Y_network(n_nodes-n_DG+i_DG,2) == 2  %it is a GFM DG
        i_GFM = i_GFM + 1;
        [AGFM0,BGFM0,CGFM0,DGFM0] = GFM_subs(Y_DER,bus_sol,GFM_Data,i_GFM,i_DG,A_GFM0,B_GFM0,C_GFM0,D_GFM0);
        A_DG_0{i_DG,i_DG}   = AGFM0;
        B_DG_0{i_DG,i_DG}   = BGFM0;
        C_DG_0{i_DG,i_DG}   = CGFM0;
        D_DG_0{i_DG,i_DG}   = DGFM0;
    elseif Y_network(n_nodes-n_DG+i_DG,2) == 3 %it is a GFL DG
        i_GFL = i_GFL + 1;
        [AGFL0,BGFL0,CGFL0,DGFL0] = GFL_subs(Y_DER,bus_sol,GFL_Data,i_GFL,i_DG,A_GFL0,B_GFL0,C_GFL0,D_GFL0);
        A_DG_0{i_DG,i_DG}   = AGFL0;
        B_DG_0{i_DG,i_DG}   = BGFL0;
        C_DG_0{i_DG,i_DG}   = CGFL0;
        D_DG_0{i_DG,i_DG}   = DGFL0;
    elseif Y_network(n_nodes-n_DG+i_DG,2) == 4  %it is a SM DG
        i_SG = i_SG + 1;
        [ASG0,BSG0,CSG0,DSG0] = SG_subs(Y_DER,bus_sol,SG_Data,i_SG,i_DG,A_SG0,B_SG0,C_SG0,D_SG0);
        A_DG_0{i_DG,i_DG}   = ASG0;
        B_DG_0{i_DG,i_DG}   = BSG0;
        C_DG_0{i_DG,i_DG}   = CSG0;
        D_DG_0{i_DG,i_DG}   = DSG0;
    else 
        A_DG_0{i_DG,i_DG}   = [];
        B_DG_0{i_DG,i_DG}   = [];
        C_DG_0{i_DG,i_DG}   = [];
        D_DG_0{i_DG,i_DG}   = [];
    end
end

%%%%%%%% Substituting in the nodes %%%%%%%%%
init_states_nodes = zeros(1,2*(n_nodes-n_DG));
col_index = 0;
A_Node_0 = [];
B_Node_0 = [];
C_Node_0 = [];
D_Node_0 = [];
for i_nodes = 1:(n_nodes - n_DG)
    
    wb = Parameters.wb;
    Cl = Y_line(1,6);
    
    % initializations:
    wg_0 = 1;
    vgj_lf = bus_sol(i_nodes,2);
    dj_lf = bus_sol(i_nodes,3)*pi/180;
    vgd_g_0 = real(vgj_lf*exp(1i*(dj_lf-theta_g_0)));
    vgq_g_0 = imag(vgj_lf*exp(1i*(dj_lf-theta_g_0)));
    
    init_states_nodes(1,col_index+1:col_index+2) = [vgd_g_0; vgq_g_0];
    col_index = col_index + 2;
    
    ishd_g_0   = -wg_0*Cl*vgq_g_0;
    ishq_g_0   =  wg_0*Cl*vgd_g_0;
    
    A_Node_0{i_nodes}   = double(vpa(subs(A_Node0)));
    B_Node_0{i_nodes}   = double(vpa(subs(B_Node0)));
    C_Node_0{i_nodes}   = double(vpa(subs(C_Node0)));
    D_Node_0{i_nodes}   = double(vpa(subs(D_Node0)));
end

%%%%%%%% Substituting in the lines %%%%%%%%%
init_states_lines = zeros(1,2*(n_lines));
col_index = 0;
A_Line_0 = [];
B_Line_0 = [];
C_Line_0 = [];
D_Line_0 = [];
for i_lines = 1:n_lines
    
    wb = Parameters.wb;
    Rl = Y_line(i_lines,4);
    Ll = Y_line(i_lines,5);
    
    % initializations:
    wg_0 = 1;
    vj_lf = bus_sol(Y_line(i_lines,2),2);
    vk_lf = bus_sol(Y_line(i_lines,3),2);
    dj_lf = bus_sol(Y_line(i_lines,2),3)*pi/180;
    dk_lf = bus_sol(Y_line(i_lines,3),3)*pi/180;
    vgdj_g_0 = real(vj_lf*exp(1i*(dj_lf-theta_g_0)));
    vgqj_g_0 = imag(vj_lf*exp(1i*(dj_lf-theta_g_0)));
    vgdk_g_0 = real(vk_lf*exp(1i*(dk_lf-theta_g_0)));
    vgqk_g_0 = imag(vk_lf*exp(1i*(dk_lf-theta_g_0)));
    
    syms ild_g_0 ilq_g_0
    E1       = [vgdj_g_0 - vgdk_g_0 - Rl*ild_g_0 + wg_0*Ll*ilq_g_0 ==0, vgqj_g_0 - vgqk_g_0 - Rl*ilq_g_0 - wg_0*Ll*ild_g_0 ==0];
    S1       = solve(E1,ild_g_0,ilq_g_0);
    ild_g_0    = double(S1.ild_g_0);
    ilq_g_0    = double(S1.ilq_g_0);
    
    init_states_lines(1,col_index+1:col_index+2) = [ild_g_0; ilq_g_0];
    col_index = col_index + 2;
    
    A_Line_0{i_lines}   = double(vpa(subs(A_Line0)));
    B_Line_0{i_lines}   = double(vpa(subs(B_Line0)));
    C_Line_0{i_lines}   = double(vpa(subs(C_Line0)));
    D_Line_0{i_lines}   = double(vpa(subs(D_Line0)));
end

%%%%%%%% Substituting in the loads %%%%%%%%%
lds = 0;
init_states_loads = zeros(1,2*(n_loads));
col_index = 0;
A_Load_0 = [];
B_Load_0 = [];
C_Load_0 = [];
D_Load_0 = [];
for i_loads = 1:n_nodes
    
    wb = Parameters.wb;
    if Y_network(i_loads,3) ~= 0
        lds = lds + 1;
        Rc = Load.R_pu(lds);
        Lc = Load.L_pu(lds);
    else
        continue
    end

    % initializations:
    wg_0 = 1;
    vgj_lf = bus_sol(i_loads,2);
    dj_lf = bus_sol(i_loads,3)*pi/180;
    vgd_g_0 = real(vgj_lf*exp(1i*(dj_lf-theta_g_0)));
    vgq_g_0 = imag(vgj_lf*exp(1i*(dj_lf-theta_g_0)));
    
    syms icd_g_0 icq_g_0
    E1       = [vgd_g_0 - Rc*icd_g_0 + wg_0*Lc*icq_g_0 ==0, vgq_g_0 - Rc*icq_g_0 - wg_0*Lc*icd_g_0 ==0];
    S1       = solve(E1,icd_g_0,icq_g_0);
    icd_g_0    = double(S1.icd_g_0);
    icq_g_0    = double(S1.icq_g_0);
    
    init_states_loads(1,col_index+1:col_index+2) = [icd_g_0; icq_g_0];
    col_index = col_index + 2;
    
    A_Load_0{i_loads}   = double(vpa(subs(A_Load0)));
    B_Load_0{i_loads}   = double(vpa(subs(B_Load0)));
    C_Load_0{i_loads}   = double(vpa(subs(C_Load0)));
    D_Load_0{i_loads}   = double(vpa(subs(D_Load0)));
end
% Find rows with zeros
rows_with_zeros = any(init_states_loads == 0, 2);
% Remove rows with zeros
init_states_loads = init_states_loads(~rows_with_zeros, :);

%% Concatenating the Individual State Spaces:

A_ol      = blkdiag (A_s_0{:},A_DG_0{:},A_Node_0{:},A_Line_0{:},A_Load_0{:});
B_ol      = blkdiag (B_s_0{:},B_DG_0{:},B_Node_0{:},B_Line_0{:},B_Load_0{:});
C_ol      = blkdiag (C_s_0{:},C_DG_0{:},C_Node_0{:},C_Line_0{:},C_Load_0{:});
D_ol      = blkdiag (D_s_0{:},D_DG_0{:},D_Node_0{:},D_Line_0{:},D_Load_0{:});

n_total_states    = (sflag_IB*nStates_s + sflag_GFM*nStates_s + sflag_SG*nStates_s) + (n_GFM - sflag_GFM)*nStates_GFM + (n_SG - sflag_SG)*nStates_SG + n_GFL *nStates_GFL + (n_nodes-n_DG)*nStates_Node + n_loads*nStates_Load + n_lines*nStates_Line;
n_total_inputs    = (sflag_IB*nIP_s(1) + sflag_GFM*nIP_s(1) + sflag_SG*nIP_s(1)) + (n_GFM - sflag_GFM)*nIP_GFM(1) + (n_SG - sflag_SG)*nIP_SG(1) + n_GFL *nIP_GFL(1) + (n_nodes-n_DG)*nIP_Node(1) + n_loads*nIP_Load(1) + n_lines*nIP_Line(1);
n_total_inputs_s  = (sflag_IB*nIP_s(2) + sflag_GFM*nIP_s(2) + sflag_SG*nIP_s(2)) + (n_GFM - sflag_GFM)*nIP_GFM(2) + (n_SG - sflag_SG)*nIP_SG(2) + n_GFL *nIP_GFL(2) + (n_nodes-n_DG)*nIP_Node(2) + n_loads*nIP_Load(2) + n_lines*nIP_Line(2);
n_total_inputs_g  = (sflag_IB*nIP_s(3) + sflag_GFM*nIP_s(3) + sflag_SG*nIP_s(3)) + (n_GFM - sflag_GFM)*nIP_GFM(3) + (n_SG - sflag_SG)*nIP_SG(3) + n_GFL *nIP_GFL(3) + (n_nodes-n_DG)*nIP_Node(3) + n_loads*nIP_Load(3) + n_lines*nIP_Line(3);
n_total_outputs   = (sflag_IB*nOP_s(1) + sflag_GFM*nOP_s(1) + sflag_SG*nOP_s(1)) + (n_GFM - sflag_GFM)*nOP_GFM(1) + (n_SG - sflag_SG)*nOP_SG(1) + n_GFL *nOP_GFL(1) + (n_nodes-n_DG)*nOP_Node(1) + n_loads*nOP_Load(1) + n_lines*nOP_Line(1);
n_total_outputs_s = (sflag_IB*nOP_s(2) + sflag_GFM*nOP_s(2) + sflag_SG*nOP_s(2)) + (n_GFM - sflag_GFM)*nOP_GFM(2) + (n_SG - sflag_SG)*nOP_SG(2) + n_GFL *nOP_GFL(2) + (n_nodes-n_DG)*nOP_Node(2) + n_loads*nOP_Load(2) + n_lines*nOP_Line(2);
n_total_outputs_g = (sflag_IB*nOP_s(3) + sflag_GFM*nOP_s(3) + sflag_SG*nOP_s(3)) + (n_GFM - sflag_GFM)*nOP_GFM(3) + (n_SG - sflag_SG)*nOP_SG(3) + n_GFL *nOP_GFL(3) + (n_nodes-n_DG)*nOP_Node(3) + n_loads*nOP_Load(3) + n_lines*nOP_Line(3);

[F,G,K,L] = assoc_matrices (N_s,N_G,N_CF,N_Cf,Y_network,Y_line1,n_total_inputs,n_total_inputs_s,n_total_inputs_g,n_total_outputs,n_total_outputs_s,n_total_outputs_g,n_nodes,n_DG,n_GFM,n_GFL,n_SG,nIP_s,nIP_GFL,nIP_GFM,nIP_SG,nIP_Node,nIP_Line,nIP_Load,sflag_GFM,sflag_SG,nOP_s,nOP_GFM,nOP_GFL,nOP_SG,nOP_Line,nOP_Load,nOP_Node,n_lines,n_loads);

n_EYE     = length((D_ol*G));
E_ol      = inv(eye(n_EYE,n_EYE) - D_ol*G);
 
A_tot0    = A_ol + B_ol*G*E_ol*C_ol;
B_tot0    = (B_ol*G*E_ol*D_ol + B_ol)*F;
C_tot0    = L*E_ol*C_ol;
D_tot0    = L*E_ol*D_ol*F + K;

[inputNames,outputNames,stateNames] = ios_names(n_nodes,n_DG,inputs_s,inputs_GFL,inputs_GFM,inputs_SG,outputs_s,outputs_GFL,outputs_GFM,outputs_SG,outputs_Node,outputs_Line,outputs_Load,n_total_inputs_s,n_total_outputs_s,n_GFM,n_GFL,n_SG,nIP_GFL,nIP_GFM,nIP_SG,nIP_Node,sflag_GFM,sflag_SG,nOP_GFM,nOP_GFL,nOP_SG,nOP_Line,nOP_Load,nOP_Node,states_s,states_GFL,states_GFM,states_SG,n_lines,n_loads,states_Line,states_Load,states_Node,nStates_s,nStates_SG,nStates_GFM,nStates_GFL,nStates_Line,nStates_Load,nStates_Node);

%% EigenValue Calculation:
[Vu,Du,Wu] = eig(A_tot0,Balance_option);
% Sorted order:
[d,ind] = sort(diag(Du));
eigA = d; %eigenVal are the sorted eigenvalues 
unstableEig = find(real(eigA) > 0);
eigA(unstableEig)

%% Modal Analysis:  
% Modal Analysis Data inputs:
% General: modal_analysis(plot_order,A_tot0,B_tot0,C_tot0,D_tot0,stateNames,inputNames,outputNames,mode_sensitivity,mode_participation,mode_mshape,n_exc_ln_ld,start_from_eig,end_eig,StepResp_Data,FreeResp_Data,init_V,xz,yz)
% plot_order = [sens_plot, part_matrix_plot, single_mode_part, eig_plot, mode_shape, free_motion, step_resp] 
% StepResp_Data = [step_offset step_amplitude Tf stepT input output] 
% FreeResp_Data = [Tf timeStep free_resp_state]

latexTab = 0;
selected_modes = [19,31,33,35,38]; %max 5 modes
% zoomed in values for eig_val plot:
xz = [-300 5];
yz = [-55 55];

n_exc_ln_ld = nStates_s+nStates_SG*(n_SG-sflag_SG)+nStates_GFM*(n_GFM-sflag_GFM)+nStates_GFL*n_GFL;
start_from_eig = 41;
end_eig = 78;
n_exc_ln_ld = n_total_states;
start_from_eig = 1;
end_eig = length(A_tot0);
StepResp_Data = [0 StepM 5 1e-3 input output];
FreeResp_Data = [5 1e-3 free_resp_state];
plot_order = [sensitivity,Participation,single_mode,EigenVal_Plot,mode_shape,Free_motion,step_resp];

% concatenate all init_states:
init_V = zeros(1,n_total_states);
if sflag_IB ==1 
    init_V(1,1:nStates_s) = init_states_IB;
elseif sflag_GFM ==1
    init_V(1,1:nStates_s) = init_states_GFM(1,:);
elseif sflag_SG ==1
    init_V(1,1:nStates_s) = init_states_SG(1,:);
end
GFm = sflag_GFM; GFl = 0; Sg = sflag_SG;
col_index = nStates_s;
for i =1:n_DG-1
    if Y_network(i+1,2) == 2 % a GFM
        GFm = GFm + 1;
        init_V(1,col_index+1:col_index+nStates_GFM) = init_states_GFM(GFm,:);
        col_index = col_index +nStates_GFM; 
    elseif Y_network(i+1,2) == 3 % a GFL
        GFl = GFl + 1;
        init_V(1,col_index+1:col_index+nStates_GFL) = init_states_GFL(GFl,:);
        col_index = col_index +nStates_GFL;        
    elseif Y_network(i+1,2) == 4 % a SG
        Sg = Sg + 1;
        init_V(1,col_index+1:col_index+nStates_SG) = init_states_SG(Sg,:);
        col_index = col_index +nStates_SG;        
    end
end
col_index = nStates_s+nStates_SG*(n_SG-sflag_SG)+nStates_GFM*(n_GFM-sflag_GFM)+nStates_GFL*n_GFL;
init_V(1,col_index+1:col_index+(n_nodes-n_DG)*nStates_Node + n_loads*nStates_Load + n_lines*nStates_Line) = [init_states_nodes init_states_lines init_states_loads];

% Finally plotting all required figures:
modal_analysis(plot_order,A_tot0,B_tot0,C_tot0,D_tot0,stateNames,inputNames,outputNames,mode_sensitivity,mode_participation,mode_mshape,n_exc_ln_ld,start_from_eig,end_eig,StepResp_Data,FreeResp_Data,init_V,xz,yz,Balance_option,latexTab,selected_modes);

end