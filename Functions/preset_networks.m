function [] = preset_networks(net_name)
% there is the following preset networks for tests:
% SM_SMIB:
% GFM_SMIB:
% GFL_SMIB:
% WSCC_3SM:
% WSCC_1SM_1GFM_1GFL:
% CIGRE_Interconnected_1SM_1GFM_1GFL:
% CIGRE_Islanded_1SM_1GFM_1GFL:

switch net_name
    case 'SM_SMIB'
    case 'GFM_SMIB'    
    case 'GFL_SMIB'
    case 'WSCC_3SM'    
    case 'WSCC_1SM_1GFM_1GFL'    
    case 'CIGRE_Interconnected_1SM_1GFM_1GFL'    
    case 'CIGRE_Islanded_1SM_1GFM_1GFL' 
        % Only for CIGRE:
        switches = 'All_on';  % Switches: All_off, All_on, S1_only, S2_only_S3_only
        % Switches:
        Tsw_S1 = 50;
        Tsw_S2 = 50;
        Tsw_S3 = 50;
        %% ========================================================================
        %  Network Description & Load Flow:
        %  ========================================================================
        % General network parameters:
        Parameters.fb           = 50; %in Hz
        Parameters.wb           = 2*pi*Parameters.fb; 
        Parameters.Sb           = 5e6; % in VA
        Parameters.Un_LV        = 10e3; % in V
        Parameters.Un_MV        = 20e3; % in V
        Parameters.Un_HV        = 110e3; % in V
        %----------------------------------------------------------%
        % Per unit calculations:
        % Since we have two voltage levels we will have two base values (one for
        % low voltage level and the other for high voltage level):
        [Ub_LV,Ib_LV,Zb_LV,Lb_LV,Cb_LV,Ub_LV_dc,Ib_LV_dc,Zb_LV_dc,Cb_LV_dc] = PU_calc(Parameters.wb,Parameters.Un_LV,Parameters.Sb);
        [Ub_MV,Ib_MV,Zb_MV,Lb_MV,Cb_MV,Ub_MV_dc,Ib_MV_dc,Zb_MV_dc,Cb_MV_dc] = PU_calc(Parameters.wb,Parameters.Un_MV,Parameters.Sb);
        [Ub_HV,Ib_HV,Zb_HV,Lb_HV,Cb_HV,Ub_HV_dc,Ib_HV_dc,Zb_HV_dc,Cb_HV_dc] = PU_calc(Parameters.wb,Parameters.Un_HV,Parameters.Sb);
        %----------------------------------------------------------%
        % Now we need to express the lines and loads in respective PU base (use LV
        % base in case the line/load is connected to LV and HV otherwise):
        % Lines Parameters:
        Line.R                 = [0.501 0.501 0.501 0.501 0.501 0.501 0.501 0.501 0.501 0.501 0.51 0.51 0.51 0.501 0.501 0.016 0.016 0.016]; 
        Line.XL                = [0.716 0.716 0.716 0.716 0.716 0.716 0.716 0.716 0.716 0.716 0.366 0.366 0.366 0.716 0.716 1.92 1.92 1.92]; 
        Line.B                 = [47.493 47.493 47.493 47.493 47.493 47.493 47.493 47.493 47.493 47.493 3.172 3.172 3.172 47.493 47.493 0 0 0].*1e-6;
        Line.C                 = Line.B ./Parameters.wb;
        Line.length            = [2.82 4.42 0.61 0.56 1.54 1.3 1.67 0.24 0.77 0.32 4.89 2.99 2 0.24 0.49 1 1 1]; %km
        Line.L                 = Line.XL ./Parameters.wb; %H/km
        Line.R_pu              = (Line.R .* Line.length) ./ Zb_MV;
        Line.L_pu              = (Line.L .* Line.length) ./ Lb_MV;
        Line.C_pu              = (Line.B .* Line.length).* Zb_MV; 
        % Loads Parameters:
        Load.P                 = [1.2     0.55      0.445    0.75     0.565   0.605    0.09     0.675    0.57     0.34     0.8     0.04     0.605].*1e6;
        Load.Q                 = [0.35    0.26638   0.11153  0.18797  0.1416  0.15163  0.05578  0.41833  0.27606  0.08521  0.2050  0.02479  0.29301].*1e6;
        Load.P_pu              = Load.P ./ Parameters.Sb;
        Load.Q_pu              = Load.Q ./ Parameters.Sb;

        %----------------------------------------------------------_
        % Load Flow:
        % GFM_controls = ["Droop","dVOC","Matching","VSM"];
        % for i_test1 = 1:length(GFM_controls)
        % DG-type: 0=none, 1=IB, 2=GFM, 3=GFL, 4=SM
        % Load = 0 -> no load connected to this node , load = 1 -> a load is
        % connected to this node
        % BusType: 1=slack bus, 2=PV bus , 3=PQ bus
        % subType: For GFM: "Droop","Droop+PLL","Droop+VIM", "dVOC",
        % "Matching", "VSM". For GFL: "PLL", "VIM". For SM/IB leave blank "". 
        % DCType: "Ideal" (GFM or GFL), "BESS" (in case of GFM),or "PV" (in case of GFL)
        %            Node   DG-type  Load  BusType  V_init  delta_init   P_gen   Q_gen   P_cons         Q_cons 
        Y_network = [1      4        0     1        1       0            0       0       0.02           0;
                     2      2        0     2        1       0            0.5     0       0              0;
                     3      3        0     3        1       0            0.4     0       0              0;
                     4      0        1     3        1       0            0       0       Load.P_pu(1)   Load.Q_pu(1);
                     5      0        1     3        1       0            0       0       Load.P_pu(2)   Load.Q_pu(2);
                     6      0        1     3        1       0            0       0       Load.P_pu(3)   Load.Q_pu(3);
                     7      0        1     3        1       0            0       0       Load.P_pu(4)   Load.Q_pu(4);
                     8      0        1     3        1       0            0       0       Load.P_pu(5)   Load.Q_pu(5);
                     9      0        1     3        1       0            0       0       Load.P_pu(6)   Load.Q_pu(6);
                     10     0        1     3        1       0            0       0       Load.P_pu(7)   Load.Q_pu(7);
                     11     0        1     3        1       0            0       0       Load.P_pu(8)   Load.Q_pu(8);
                     12     0        1     3        1       0            0       0       Load.P_pu(9)   Load.Q_pu(9);
                     13     0        1     3        1       0            0       0       Load.P_pu(10)  Load.Q_pu(10);
                     14     0        1     3        1       0            0       0       Load.P_pu(11)  Load.Q_pu(11);
                     15     0        0     3        1       0            0       0       0              0;
                     16     0        1     3        1       0            0       0       Load.P_pu(12)  Load.Q_pu(12);
                     17     0        1     3        1       0            0       0       Load.P_pu(13)  Load.Q_pu(13)];

        % Line Description:
        %            Line  From To  R(pu)         L(pu)         C(pu)           Length(km)
        Y_line    = [1     15   14  Line.R_pu(1)  Line.L_pu(1)  Line.C_pu(1)    Line.length(1);
                     2     15   4   Line.R_pu(2)  Line.L_pu(2)  Line.C_pu(2)    Line.length(2);
                     3     7    4   Line.R_pu(3)  Line.L_pu(3)  Line.C_pu(3)    Line.length(3);
                     4     7    8   Line.R_pu(4)  Line.L_pu(4)  Line.C_pu(4)    Line.length(4);
                     5     9    8   Line.R_pu(5)  Line.L_pu(5)  Line.C_pu(5)    Line.length(5);
                     6     11   4   Line.R_pu(6)  Line.L_pu(6)  Line.C_pu(6)    Line.length(6);
                     7     11   10  Line.R_pu(7)  Line.L_pu(7)  Line.C_pu(7)    Line.length(7);
                     8     11   6   Line.R_pu(8)  Line.L_pu(8)  Line.C_pu(8)    Line.length(8);
                     9     12   6   Line.R_pu(9)  Line.L_pu(9)  Line.C_pu(9)    Line.length(9);
                     10    12   13  Line.R_pu(10) Line.L_pu(10) Line.C_pu(10)   Line.length(10);
                     11    17   16  Line.R_pu(11) Line.L_pu(11) Line.C_pu(11)   Line.length(11);
                     12    17   5   Line.R_pu(12) Line.L_pu(12) Line.C_pu(12)   Line.length(12);
                     13    11   5   Line.R_pu(13) Line.L_pu(13) Line.C_pu(13)   Line.length(13); %S1
                     14    9    10  Line.R_pu(14) Line.L_pu(14) Line.C_pu(14)   Line.length(14); %S2
                     15    7    13  Line.R_pu(15) Line.L_pu(15) Line.C_pu(15)   Line.length(15)]; %S3

        % Y_DG is complementary to Y_network:
        %            Node DG-type  subType    DCType   
        Y_DG      = [1    4        "SM"       "NA";
                     2    2        "Droop"    "Ideal";
                     3    3        "PLL"      "Ideal";
                     4    0        "NA"       "NA";
                     5    0        "NA"       "NA";
                     6    0        "NA"       "NA";
                     7    0        "NA"       "NA";
                     8    0        "NA"       "NA";
                     9    0        "NA"       "NA";
                     10   0        "NA"       "NA";
                     11   0        "NA"       "NA";
                     12   0        "NA"       "NA";
                     13   0        "NA"       "NA";
                     14   0        "NA"       "NA";
                     15   0        "NA"       "NA";
                     16   0        "NA"       "NA";
                     17   0        "NA"       "NA"];

        % Transformer Description:
        %            TR    From To    R(pu)         L(pu)         C(pu)           put "1"
        Y_TR      = [1     1    4     Line.R_pu(16) Line.L_pu(16) Line.C_pu(16)   Line.length(16);
                     2     2    5     Line.R_pu(17) Line.L_pu(17) Line.C_pu(17)   Line.length(17);
                     3     3    6     Line.R_pu(18) Line.L_pu(18) Line.C_pu(18)   Line.length(18)];
                 
        % Selection of Different Elements Voltage Level & Nominal Power:
        Nom_Power     = [Parameters.Sb; Parameters.Sb; Parameters.Sb; Parameters.Sb; Parameters.Sb; Parameters.Sb]; %1: IB (if any), 2: SMs, 3: TRs, 4: Lines, 5: Loads
        Voltage_level = ["HV"; "LV"; "MV"; "LV"; "MV"; "MV"; "LV"]; %1: IB (if any), 2: SMs, 3: TRs (prim. side), 4: TRs (sec. side), 5: Lines, 6: Loads

                 
end

        
        
end