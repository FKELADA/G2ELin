function [bus_sol,P_loss,I_line,iter,conv_flag,Y_bus] = Load_flow(Y_network, Y_line,Sb_MW)
%%%%%%%%%%%%%%%%%%%%
n_buses = size(Y_network,1);
n_lines = size(Y_line,1);
%Load Flow:
% Type: (1 for slack, 2 for PV and 3 for PQ)
%       no              V_init          delta_init      P_gen           Q_gen            P_cons           Q_cons           -                 -               Type
Bus = [ Y_network(:,1)  Y_network(:,5)  Y_network(:,6)  Y_network(:,7)  Y_network(:,8)  Y_network(:,9)  Y_network(:,10)  zeros(n_buses,1)  zeros(n_buses,1) Y_network(:,4)]; 

%       From          To            R_line(pu)    L_line(pu)     B_line(pu)
Line = [Y_line(:,2)   Y_line(:,3)   Y_line(:,4)   Y_line(:,5)    Y_line(:,6)];   

[bus_sol,P_loss,I_line,iter,conv_flag,Y_bus] = Power_Fl(Line,Bus,Sb_MW);


% data = [bus_sol(:,1)  bus_sol(:,10)  bus_sol(:,2)  bus_sol(:,3).*(pi/180)  bus_sol(:,4)  bus_sol(:,5)  bus_sol(:,6)  bus_sol(:,7)] ;
% CN   = {'Bus'  'Type'  'Voltage(pu)'  'Angle(rad)'  'P_gen(pu)'  'Q_gen(pu)'  'P_cons(pu)'  'Q_cons(pu)'};
% figure , t = uitable('Data', data, 'ColumnName', CN ,'Position', [25.599999999999994,210.6,523.4,79.79999999999973]);
% 



% % Y Bus Formation
% a = size(Bus); n1 = a(1);
% b = size(Line) ; n2 = b(1) ;
% Y_Bus=zeros(n1);
% for I = 1:n1
%     for J = 1:n1
%         if I == J
%             for m = 1 : n2
%                 if Line(m,1)==I | Line(m,2)==I
%                     Y_Bus(I,I) = Y_Bus(I,I)+ (1/(Line(m,3)+Line(m,4)*1i))+Line(m,5);
%                 end
%             end
%         else
%             for m = 1 : n2
%                 if (Line(m,1)==I & Line(m,2)==J ) | (Line(m,1)==J & Line(m,2)==I )
%                     Y_Bus(I,J) = Y_Bus(I,J)+ (-1/((Line(m,3)+Line(m,4)*1i)));break
%                 end
%             end
%         end
%     end
% end
% 
% Y_Bus
% 
% P_nodes = zeros(n_buses,1);
% V_vector  = zeros(n_buses,1);
% for I  = 1:n_buses
%     V_vector(I) = bus_sol(I,2)*cos(bus_sol(I,3).*(pi/180)) + 1i*bus_sol(I,2)*sin(bus_sol(I,3).*(pi/180));
% end
% for I = 1:n_buses
%     [row,col] = find(Y_line(:,2:3)==I);
%     for J = 1:numel(row)
%         
%     end
%     P_nodes(I) = real(V_vector(Y_line(I,2))*conj(Y_Bus(Y_line(I,2),Y_line(I,3)))*V_vector(Y_line(I,3)));
% end

end