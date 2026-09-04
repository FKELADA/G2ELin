function [results_matrix] = SCR2(Y_bus,bus_sol,Y_DER)
c = 0 ;% counter for the iteration number.
c1 = 1 ;% counter to help find parallel bus
n_linhas = size(Y_line,1);
bus = bus_sol(:,1) ;  % bus
V = bus_sol(:,2) ;    % voltage
theta = bus_sol(:,3) ;% angle in degrees
type = bus_sol(:,10) ; % 1 se Vtheta , 2 se PV , 3 se PQ
Pg = bus_sol(:,4) ;   % Active power generated in the bus
Qg = bus_sol(:,5) ;   % Reactive power generated in the bus
Xd = Y_DER(:,13) ;   % Generator reactance
Pl = bus_sol(:,6) ;   % Active power delivered to the bus
Ql = bus_sol(:,7) ;   % Reactive power delivered to the bus
Sl = Pl + Ql*1j ;      % Complex power delivered to the bus
f1 = find(type == 1) ; Vtheta = f1' ; % angle reference bus
f2 = find(type == 2) ; PV = f2' ; % buses with load and generation or synchronous compensator
f3 = find(type == 3) ; PQ = f3' ; % buses with load
%% Memory Pre-allocation
nbus = length(bus) ; % number of system buses
nDER = size(Y_DER,1);
Yd = zeros(nbus) ;           % calculated admittance with generator reactance
Zl = zeros(1,nbus) ;         % load impedance obtained by Z=V²/S*
If = zeros(nbus,1) ;         % column vector of fault currents
cc = 1 ;
for c = 1 : nDER
     V(c) = V(c)*((cosd(theta(c))+1j*sind(theta(c)))) ;
    if Xd(c) ~= 0
       Yd(c,c) = 1/(Xd(c)*1j) ; % admittance obtained by the generator reactance
    end
end
%%                  Building the Y_bus Matrix
Y_bus = Y_bus + Yd ; % total system admittance
G = real(Y_bus) ;   % conductance is the real part of the admittance matrix.
B = imag(Y_bus) ;   % susceptance is the imaginary part of the admittance matrix
Zbus = inv(Y_bus) ; 
If(cc) = - V(cc)/Zbus(cc,cc) ;
deltaV = Y_bus \ If ; % Zbus * If
Vf = V + deltaV ;  % voltage V during fault
Vabs = abs(Vf) ;   % phase voltage module
thetaV = angle(Vf) * 180 /pi ;
If_fault = zeros(nbus) ; 
for il = 1 : n_linhas    % ranging from 1 to the number of lines.
    k = Y_line(il,2) ;% origin bus.
    m = Y_line(il,3) ;% destination bus.
    Rlin = Y_line(il,4) ;% resistance 
    Xlin = (Y_line(il,5)*1j) ;% reactance
    z = Rlin + Xlin ;
    If_fault(k,m) = (Vf(k) - Vf(m)) / z ; 
end

end