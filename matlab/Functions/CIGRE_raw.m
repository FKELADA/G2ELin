function [Y_network1,Y_line1,Y_TR1,nodes_from_r,nodes_to_r,Line_RX,Load_P,Load_Q,Line_R,Line_L,Line_C,Line_length] = CIGRE_raw(S0,S1,S2,S3,wb,Sb,Zb_lines,Lb_lines)
% This script describes the raw form of the CIGRE MV benchmark network
% S0: if 1 means we are using the interconnected mode, else it is the
% islanded mode
% S1: it is the switch between node 8 & node 12
% S2: it is the switch between node 6 & node 7
% S3: it is the switch between node 4 & node 11

% Lines Parameters:
Line.R                 = [0.501   0.501   0.501   0.501   0.501   0.501   0.501   0.501   0.501   0.501   0.51   0.51   0.51   0.501   0.501]; 
Line.XL                = [0.716   0.716   0.716   0.716   0.716   0.716   0.716   0.716   0.716   0.716   0.366  0.366  0.366  0.716   0.716]; 
Line.B                 = [47.493  47.493  47.493  47.493  47.493  47.493  47.493  47.493  47.493  47.493  3.172  3.172  3.172  47.493  47.493].*1e-6;
Line.length            = [2.82    4.42    0.61    0.56    1.54    1.3     1.67    0.24    0.77    0.32    4.89   2.99   2      0.24    0.49]; %km
Line.C                 = Line.B ./wb;
Line.L                 = Line.XL ./wb; %H/km
Line.R_pu              = (Line.R .* Line.length) ./ Zb_lines;
Line.L_pu              = (Line.L .* Line.length) ./ Lb_lines;
Line.C_pu              = (Line.B .* Line.length).* Zb_lines;
Line_RX                = Line.R_pu ./ Line.L_pu;

% Loads Parameters:
Load.P                 = [1.2    0     0.55      0.445     0.34      0.57      0.675     0.605      0.09       0.565     0.75      0.605    0.04     0.8   ].*1e6;
Load.Q                 = [0.35   0     0.26638   0.11153   0.08521   0.27606   0.41833   0.15163    0.05578    0.1416    0.18797   0.29301  0.02479  0.2050].*1e6;
Load.P_pu              = Load.P ./ Sb;
Load.Q_pu              = Load.Q ./ Sb;

% Outputs:
Load_P = Load.P;
Load_Q = Load.Q;
Line_R = Line.R;
Line_L = Line.L;
Line_C = Line.C;
Line_length = Line.length;

if S0 == 0  % islanded mode
% GFM_controls = ["Droop","dVOC","Matching","VSM"];
% for i_test1 = 1:length(GFM_controls)
% DG-type: 0=none, 1=IB, 2=GFM, 3=GFL, 4=SM
% Load = 0 -> no load connected to this node , load = 1 -> a load is
% connected to this node
% BusType: 1=slack bus, 2=PV bus , 3=PQ bus
% subType: For GFM: "Droop","Droop+PLL","Droop+VIM", "dVOC",
% "Matching", "VSM". For GFL: "PLL", "VIM". For SM/IB leave blank "". 
% DCType: "Ideal" (GFM or GFL), "BESS" (in case of GFM),or "PV" (in case of GFL)
%            Node    DG-type  Load  BusType  V_init  delta_init   P_gen   Q_gen   P_cons         Q_cons 
Y_network1 = [1      0        1     3        1       0            0       0       Load.P_pu(1)   Load.Q_pu(1);
              2      0        0     3        1       0            0       0       Load.P_pu(2)   Load.Q_pu(2);
              3      0        1     3        1       0            0       0       Load.P_pu(3)   Load.Q_pu(3);
              4      0        1     3        1       0            0       0       Load.P_pu(4)   Load.Q_pu(4);
              5      0        1     3        1       0            0       0       Load.P_pu(5)   Load.Q_pu(5);
              6      0        1     3        1       0            0       0       Load.P_pu(6)   Load.Q_pu(6);
              7      0        1     3        1       0            0       0       Load.P_pu(7)   Load.Q_pu(7);
              8      0        1     3        1       0            0       0       Load.P_pu(8)   Load.Q_pu(8);
              9      0        1     3        1       0            0       0       Load.P_pu(9)   Load.Q_pu(9);
              10     0        1     3        1       0            0       0       Load.P_pu(10)  Load.Q_pu(10);
              11     0        1     3        1       0            0       0       Load.P_pu(11)  Load.Q_pu(11);
              12     0        1     3        1       0            0       0       Load.P_pu(12)  Load.Q_pu(12);
              13     0        1     3        1       0            0       0       Load.P_pu(13)  Load.Q_pu(13);
              14     0        1     3        1       0            0       0       Load.P_pu(14)  Load.Q_pu(14)];

          
% Line Description:
%            Line   From To  R(pu)         L(pu)         C(pu)           Length(km)
Y_line1    = [1     2    1   Line.R_pu(1)  Line.L_pu(1)  Line.C_pu(1)    Line.length(1);
              2     2    3   Line.R_pu(2)  Line.L_pu(2)  Line.C_pu(2)    Line.length(2);
              3     4    3   Line.R_pu(3)  Line.L_pu(3)  Line.C_pu(3)    Line.length(3);
              4     4    5   Line.R_pu(4)  Line.L_pu(4)  Line.C_pu(4)    Line.length(4);
              5     6    5   Line.R_pu(5)  Line.L_pu(5)  Line.C_pu(5)    Line.length(5);
              6     8    3   Line.R_pu(6)  Line.L_pu(6)  Line.C_pu(6)    Line.length(6);
              7     8    7   Line.R_pu(7)  Line.L_pu(7)  Line.C_pu(7)    Line.length(7);
              8     8    9   Line.R_pu(8)  Line.L_pu(8)  Line.C_pu(8)    Line.length(8);
              9     10   9   Line.R_pu(9)  Line.L_pu(9)  Line.C_pu(9)    Line.length(9);
              10    10   11  Line.R_pu(10) Line.L_pu(10) Line.C_pu(10)   Line.length(10);
              11    13   14  Line.R_pu(11) Line.L_pu(11) Line.C_pu(11)   Line.length(11);
              12    13   12  Line.R_pu(12) Line.L_pu(12) Line.C_pu(12)   Line.length(12);
              13.*S1    8.*S1    12.*S1  Line.R_pu(13).*S1 Line.L_pu(13).*S1 Line.C_pu(13).*S1   Line.length(13).*S1; %S1
              14.*S2    6.*S2    7.*S2   Line.R_pu(14).*S2 Line.L_pu(14).*S2 Line.C_pu(14).*S2   Line.length(14).*S2; %S2
              15.*S3    4.*S3    11.*S3  Line.R_pu(15).*S3 Line.L_pu(15).*S3 Line.C_pu(15).*S3   Line.length(15).*S3]; %S3


else       % interconnected mode
% GFM_controls = ["Droop","dVOC","Matching","VSM"];
% for i_test1 = 1:length(GFM_controls)
% DG-type: 0=none, 1=IB, 2=GFM, 3=GFL, 4=SM
% Load = 0 -> no load connected to this node , load = 1 -> a load is
% connected to this node
% BusType: 1=slack bus, 2=PV bus , 3=PQ bus
% subType: For GFM: "Droop","Droop+PLL","Droop+VIM", "dVOC",
% "Matching", "VSM". For GFL: "PLL", "VIM". For SM/IB leave blank "". 
% DCType: "Ideal" (GFM or GFL), "BESS" (in case of GFM),or "PV" (in case of GFL)
%            Node    DG-type  Load  BusType  V_init  delta_init   P_gen   Q_gen   P_cons         Q_cons 
Y_network1 = [1      0        1     3        1       0            0       0       Load.P_pu(1)+Load.P_pu(14)   Load.Q_pu(1)+Load.Q_pu(14);
              2      0        0     3        1       0            0       0       Load.P_pu(2)   Load.Q_pu(2);
              3      0        1     3        1       0            0       0       Load.P_pu(3)   Load.Q_pu(3);
              4      0        1     3        1       0            0       0       Load.P_pu(4)   Load.Q_pu(4);
              5      0        1     3        1       0            0       0       Load.P_pu(5)   Load.Q_pu(5);
              6      0        1     3        1       0            0       0       Load.P_pu(6)   Load.Q_pu(6);
              7      0        1     3        1       0            0       0       Load.P_pu(7)   Load.Q_pu(7);
              8      0        1     3        1       0            0       0       Load.P_pu(8)   Load.Q_pu(8);
              9      0        1     3        1       0            0       0       Load.P_pu(9)   Load.Q_pu(9);
              10     0        1     3        1       0            0       0       Load.P_pu(10)  Load.Q_pu(10);
              11     0        1     3        1       0            0       0       Load.P_pu(11)  Load.Q_pu(11);
              12     0        1     3        1       0            0       0       Load.P_pu(12)  Load.Q_pu(12);
              13     0        1     3        1       0            0       0       Load.P_pu(13)  Load.Q_pu(13)];

          
% Line Description:
%            Line   From To  R(pu)         L(pu)         C(pu)           Length(km)
Y_line1    = [1     2    1   Line.R_pu(1)  Line.L_pu(1)  Line.C_pu(1)    Line.length(1);
              2     2    3   Line.R_pu(2)  Line.L_pu(2)  Line.C_pu(2)    Line.length(2);
              3     4    3   Line.R_pu(3)  Line.L_pu(3)  Line.C_pu(3)    Line.length(3);
              4     4    5   Line.R_pu(4)  Line.L_pu(4)  Line.C_pu(4)    Line.length(4);
              5     6    5   Line.R_pu(5)  Line.L_pu(5)  Line.C_pu(5)    Line.length(5);
              6     8    3   Line.R_pu(6)  Line.L_pu(6)  Line.C_pu(6)    Line.length(6);
              7     8    7   Line.R_pu(7)  Line.L_pu(7)  Line.C_pu(7)    Line.length(7);
              8     8    9   Line.R_pu(8)  Line.L_pu(8)  Line.C_pu(8)    Line.length(8);
              9     10   9   Line.R_pu(9)  Line.L_pu(9)  Line.C_pu(9)    Line.length(9);
              10    10   11  Line.R_pu(10) Line.L_pu(10) Line.C_pu(10)   Line.length(10);
              11    13   1  Line.R_pu(11) Line.L_pu(11) Line.C_pu(11)   Line.length(11);
              12    13   12  Line.R_pu(12) Line.L_pu(12) Line.C_pu(12)   Line.length(12);
              13.*S1    8.*S1    12.*S1  Line.R_pu(13).*S1 Line.L_pu(13).*S1 Line.C_pu(13).*S1   Line.length(13).*S1; %S1
              14.*S2    6.*S2    7.*S2   Line.R_pu(14).*S2 Line.L_pu(14).*S2 Line.C_pu(14).*S2   Line.length(14).*S2; %S2
              15.*S3    4.*S3    11.*S3  Line.R_pu(15).*S3 Line.L_pu(15).*S3 Line.C_pu(15).*S3   Line.length(15).*S3]; %S3


end


nonZeroRows = any(Y_line1, 2);
Y_line1 = Y_line1(nonZeroRows, :);
         
% FOR CIGRE_RAW:
% nodes_from_r = [2 4 6 8 10 13];
% nodes_to_r   = [1 3 5 7 9  11 12 14];     
nodes_from_r = sort(unique(Y_line1(:,2)));
nodes_to_r   = sort(unique(Y_line1(:,3)));

Y_TR1      = [];
end