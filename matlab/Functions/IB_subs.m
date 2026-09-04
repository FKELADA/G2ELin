function [A_IB_0,B_IB_0,C_IB_0,D_IB_0] = IB_subs(Y_DER,N_s,bus_sol,IB_Data,i_DG,A_0,B_0,C_0,D_0)

    Rup       = IB_Data(1,1);
    Lup       = IB_Data(1,2);
    theta_up_0 = IB_Data(1,3);
    wup_0      = IB_Data(1,4);
    Vup_0     = IB_Data(1,5);
    P_0       = bus_sol(N_s,4);
    Q_0       = bus_sol(N_s,5);
    wb        = IB_Data(1,6);

    vgd_0     = real((bus_sol(Y_DER(1,1),2))*exp(1i*((bus_sol(Y_DER(1,1),3))*pi/180-theta_up_0)));
    vgq_0     = imag((bus_sol(Y_DER(1,1),2))*exp(1i*((bus_sol(Y_DER(1,1),3))*pi/180-theta_up_0)));
    vgd_g_0   = vgd_0;
    vgq_g_0   = vgq_0;
    igd_g_0  = real((P_0 - 1i*Q_0)/Vup_0);
    igq_g_0  = imag((P_0 - 1i*Q_0)/Vup_0);
%     syms igd_g_0 igq_g_0
%     E1       = [Vup_0 - vgd_g_0 - Rup*igd_g_0 + wup_0*Lup*igq_g_0 ==0, vgq_g_0 - Rup*igq_g_0 - wup_0*Lup*igd_g_0 ==0];
%     S1       = solve(E1,igd_g_0,igq_g_0);
%     igd_g_0    = double(S1.igd_g_0);
%     igq_g_0    = double(S1.igq_g_0);
    

    A_IB_0   = double(vpa(subs(A_0)));
    B_IB_0   = double(vpa(subs(B_0)));
    C_IB_0   = double(vpa(subs(C_0)));
    D_IB_0   = double(vpa(subs(D_0)));
end