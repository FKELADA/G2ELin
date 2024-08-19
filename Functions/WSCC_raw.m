function [Y_network1,Y_line1,Y_TR1,nodes_from_r,nodes_to_r,Line_RX,Load_P,Load_Q,Line_R,Line_L,Line_C,Line_length] = WSCC_raw(wb,Sb,Zb_lines,Lb_lines,modif)
% This script describes the raw form of the WSCC benchmark network
% Lines Parameters:
Line.R_pu              = [0.01 0.017 0.032*1 0.0085 0.0119 0.039].*modif; 
Line.L_pu              = [0.085 0.092 0.161 0.072 0.1008 0.17]; 
Line.B_pu              = [0.088 0.079 0.153 0.0745 0.1045 0.179].*2; 
Line.length            = [1 1 1 1 1 1]; %km
Line.XL                = Line.L_pu .*wb; 
Line.R                 = (Line.R_pu ./ Line.length) .* Zb_lines; %ohm/km
Line.L                 = (Line.L_pu ./ Line.length) .* Lb_lines; %H/km
Line.B                 = (Line.B_pu ./ Line.length) ./ Zb_lines; % B = omega*C = B_pu/Zb
Line.C                 = Line.B ./ (wb); % F/km
Line_RX                = Line.R_pu ./ Line.L_pu;
% Loads Parameters:
Load.P                 = [0  0  0  125 100 90].*1e6;
Load.Q                 = [0  0  0  50  35  30].*1e6;
Load.P_pu              = Load.P ./ Sb;
Load.Q_pu              = Load.Q ./ Sb;
% Outputs:
Load_P = Load.P;
Load_Q = Load.Q;
Line_R = Line.R;
Line_L = Line.L;
Line_C = Line.C;
Line_length = Line.length;
% DG-type: 0=none, 1=IB, 2=GFM, 3=GFL, 4=SM
% Load = 0 -> no load connected to this node , load = 1 -> a load is
% connected to this node
% BusType: 1=slack bus, 2=PV bus , 3=PQ bus
%            Node    DG-type  Load  BusType  V_init  delta_init   P_gen   Q_gen   P_cons         Q_cons 
Y_network1 = [1      0        0     3        1       0            0       0       Load.P_pu(1)   Load.Q_pu(1);
              2      0        0     3        1       0            0       0       Load.P_pu(2)   Load.Q_pu(2);
              3      0        0     3        1       0            0       0       Load.P_pu(3)   Load.Q_pu(3);
              4      0        1     3        1       0            0       0       Load.P_pu(4)   Load.Q_pu(4);
              5      0        1     3        1       0            0       0       Load.P_pu(5)   Load.Q_pu(5);
              6      0        1     3        1       0            0       0       Load.P_pu(6)   Load.Q_pu(6)];

          
% Line Description:
%            Line   From To  R(pu)         L(pu)         B(pu)           Length(km)
Y_line1    = [1     1    4   Line.R_pu(1)  Line.L_pu(1)  Line.B_pu(1)    Line.length(1);
              2     1    6   Line.R_pu(2)  Line.L_pu(2)  Line.B_pu(2)    Line.length(2);
              3     2    4   Line.R_pu(3)  Line.L_pu(3)  Line.B_pu(3)    Line.length(3);
              4     2    5   Line.R_pu(4)  Line.L_pu(4)  Line.B_pu(4)    Line.length(4);
              5     3    5   Line.R_pu(5)  Line.L_pu(5)  Line.B_pu(5)    Line.length(5);
              6     3    6   Line.R_pu(6)  Line.L_pu(6)  Line.B_pu(6)    Line.length(6)]; 


nonZeroRows = any(Y_line1, 2);
Y_line1 = Y_line1(nonZeroRows, :);
             
nodes_from_r = sort(unique(Y_line1(:,2)));
nodes_to_r   = sort(unique(Y_line1(:,3)));

Y_TR1      = [];

end