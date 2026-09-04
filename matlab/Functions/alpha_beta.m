function [Vab] = alpha_beta(bus_sol)
% this function computes all node voltages in the stationary reference
% frame alpha-beta (real-imaginary)
Vr = zeros(1,size(bus_sol,1));
Vi = zeros(1,size(bus_sol,1));

for i = 1:size(bus_sol,1)
    Vr(i) = bus_sol(i,2)*cos(bus_sol(i,3)*pi/180);
    Vi(i) = bus_sol(i,2)*sin(bus_sol(i,3)*pi/180);
end

Vab = [Vr.' Vi.'];
end