clc; clear all; close all
% Parameters
syms wb Rg Lg Rf Lf Cf
% States (x)
syms isd isq igd igq ved veq
%Algebraic (z)
% Inputs (u)
syms vgd vgq vmd vmq w
%%%%%%%%%%% Forming the state/input/output/algebraic vectors: %%%%%%%%%%%

algVec           = [];
inputVec         = [vgd vgq vmd vmq w];
stateVec         = [isd isq igd igq ved veq];

% Equilibrium vectors
for i = 1:length(inputVec)
    inputVecEq{i} = sym(strcat(char(inputVec(i)),'_0'));
end
inputVecEq = cell2sym(inputVecEq);

for i = 1:length(stateVec)
    stateVecEq{i} = sym(strcat(char(stateVec(i)),'_0'));
end
stateVecEq = cell2sym(stateVecEq);

% for i = 1:length(algVec)
%     algVecEq{i} = sym(strcat(char(algVec(i)),'_0'));
% end
algVecEq = [];

%%%%%%%%%%%%%%%%%%%%%%% Differential equations %%%%%%%%%%%%%%%%%%%%%%%%%%%   
disd          = (wb/Lf)*(vmd - ved - Rf*isd + w*Lf*isq);
disq          = (wb/Lf)*(vmq - veq - Rf*isq - w*Lf*isd);
digd          = (wb/Lg)*(ved - vgd - Rg*igd + w*Lg*igq);
digq          = (wb/Lg)*(veq - vgq - Rg*igq - w*Lg*igd);
dved          = (wb/Cf)*(isd - igd + w*Cf*veq);
dveq          = (wb/Cf)*(isq - igq - w*Cf*ved);

diffeqVec = [disd disq digd digq dved dveq];    

%%%%%%%%%%%%%%%%%%%%%%% Algebraic equations %%%%%%%%%%%%%%%%%%%%%%%%%%%

algeqVec    = [];

%%%%%%%%%%%%%%%%%%%%%%% Output equations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output vector

outputeqVec      = [igd igq];

%%%%%%%%%%%%%%%%%%%%%%% Linearization %%%%%%%%%%%%%%%%%%%%%%%%%%%

Fx          = simplify(jacobian(diffeqVec, stateVec));
% Fz          = simplify(jacobian(diffeqVec, algVec));
Fu          = simplify(jacobian(diffeqVec, inputVec));

Gx          = simplify(jacobian(algeqVec, stateVec));
% Gz          = simplify(jacobian(algeqVec, algVec));
Gu          = simplify(jacobian(algeqVec, inputVec));
% invGz       = inv(Gz);

Hx          = simplify(jacobian(outputeqVec, stateVec));
% Hz          = simplify(jacobian(outputeqVec, algVec));
Hu          = simplify(jacobian(outputeqVec, inputVec));

Ai          = Fx ;
Bi          = Fu ;
Ci          = Hx ;
Di          = Hu ;

A_SG0      = subs(Ai, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
B_SG0      = subs(Bi, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
C_SG0      = subs(Ci, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);
D_SG0      = subs(Di, [stateVec inputVec algVec], [stateVecEq inputVecEq algVecEq]);

