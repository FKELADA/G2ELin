function [A_SG_0,B_SG_0,C_SG_0,D_SG_0] = SG_subs(Y_DER,bus_sol,SG_Data,i_SG,i_DG,A_0,B_0,C_0,D_0)
% SG Data:
%1)wb   2)wff(pu)   3)Ra(pu)   4)Ll(pu) 5)Cg(pu) 6)Rt(pu)  7)Lt(pu)  8)Lad(pu)  9)Laq(pu)  
%10)Lfd(pu)  11)Rfd(pu)  12)L1d(pu)  13)L1q(pu)  14)L2q(pu)  15)R1d(pu) 
%16)R1q(s)   17)R2q(s)   18)H(s) 19)KD  20)mp  
%21)TG  22)T_LP  23)T_HP  24)K_PSS  25)T1n  26)T1d   
%27)T2n  28)T2d  29)Tr  30)Ka  31)Ta  32)Ke  33)Te 34)Kfd  35)Tfd 
%36)Pt  37)Qt  38)Et  39)theta0 40)delta0  41)ed_0  42)eq_0
%43)id_0 44)iq_0 45)Psi_d0 46)Psi_q0 47)Psi_fd0
%48)Psi_1d0 49)Psi_1q0 50)Psi_2q0 51) SM.Lffd 52)SM.RL_pu 53)SM.id_0 54)SM.iq_0
%55) SM.theta_0(1)

    wb        = SG_Data(i_SG,1);
    wff       = SG_Data(i_SG,2);
    Ra        = SG_Data(i_SG,3);
    Ll        = SG_Data(i_SG,4);
    Cg        = SG_Data(i_SG,5);
    Rt        = SG_Data(i_SG,6);
    Lt        = SG_Data(i_SG,7);
    Lad       = SG_Data(i_SG,8);
    Laq       = SG_Data(i_SG,9);
    Lf1d      = SG_Data(i_SG,10);
    Lffd      = SG_Data(i_SG,51);
    Lfd       = Lffd - Lad;
    Rfd       = SG_Data(i_SG,11);
    L11d      = SG_Data(i_SG,12);
    L1d       = L11d - Lad;
    L11q      = SG_Data(i_SG,13);
    L1q       = L11q - Laq;
    L22q      = SG_Data(i_SG,14);
    L2q       = L22q - Laq;
    R1d       = SG_Data(i_SG,15);
    R1q       = SG_Data(i_SG,16);
    R2q       = SG_Data(i_SG,17);
    H         = SG_Data(i_SG,18);
    KD        = SG_Data(i_SG,19);
    mp        = SG_Data(i_SG,20)/100;
    TG        = SG_Data(i_SG,21);
    T_LP      = SG_Data(i_SG,22);
    T_HP      = SG_Data(i_SG,23);
    K_PSS     = SG_Data(i_SG,24);
    T1n       = SG_Data(i_SG,25);
    T1d       = SG_Data(i_SG,26);
    T2n       = SG_Data(i_SG,27);
    T2d       = SG_Data(i_SG,28);
    Tr        = SG_Data(i_SG,29);
    Ka        = SG_Data(i_SG,30);
    Ta        = SG_Data(i_SG,31);
    Ke        = SG_Data(i_SG,32);
    Te        = SG_Data(i_SG,33);
    Kfd       = SG_Data(i_SG,34);
    Tfd       = SG_Data(i_SG,35);
    Rg        = SG_Data(i_SG,52);

    % Initial Conditions:
    % Setpoints
    p_ref_0   = SG_Data(i_SG,36);
    w_ref_0   = 1;

    
    % Inputs and state variables:
    % States Independent to others:
    theta_0   = SG_Data(i_SG,39);
    Delta_0   = SG_Data(i_SG,40);
    wr_0      = 1;
    Dwr_0     = wr_0 - w_ref_0;
    theta_g_0   = SG_Data(i_SG,55);
    ved_0     = SG_Data(i_SG,41);
    veq_0     = SG_Data(i_SG,42);
    igd_0     = SG_Data(i_SG,43);
    igq_0     = SG_Data(i_SG,44);
    id_0      = SG_Data(i_SG,53);
    iq_0      = SG_Data(i_SG,54);
    vgd_0     = real((bus_sol(Y_DER(i_DG,1),2))*exp(1i*((bus_sol(Y_DER(i_DG,1),3))*pi/180-theta_0)));
    vgq_0     = imag((bus_sol(Y_DER(i_DG,1),2))*exp(1i*((bus_sol(Y_DER(i_DG,1),3))*pi/180-theta_0)));
    vgd_g_0   = vgd_0*cos(theta_g_0 - theta_0) + vgq_0*sin(theta_g_0 - theta_0);
    vgq_g_0   = -vgd_0*sin(theta_g_0 - theta_0) + vgq_0*cos(theta_g_0 - theta_0);

    % states that depend on others:
    
    Pm_0      = p_ref_0;
    Dw1_0     = Dwr_0;
    v1_0      = K_PSS*Dw1_0;
    v2_0      = v1_0;
    vpss_0    = v2_0;
    
    Et_0      = sqrt(ved_0^2 + veq_0^2);
    e1_0      = Et_0;
    e3_0      = 0;
    ifd_0     = (veq_0 + Ra*iq_0 + (Lad+Ll)*id_0)/Lad;
    efd_0     = Lad*ifd_0;
    v_ref_0   = efd_0/Ka + e1_0;
    e2_0      = Ka*(v_ref_0 - e1_0 - e3_0 + vpss_0);

%     efd_0     = Ke*e2_0;   

    phi_d_0   = SG_Data(i_SG,45);
    phi_q_0   = SG_Data(i_SG,46);
    phi_fd_0  = SG_Data(i_SG,47);
    phi_1d_0  = SG_Data(i_SG,48);
    phi_1q_0  = SG_Data(i_SG,49);
    phi_2q_0  = SG_Data(i_SG,50);
    

    A_SG_0   = double(vpa(subs(A_0)));
    B_SG_0   = double(vpa(subs(B_0)));
    C_SG_0   = double(vpa(subs(C_0)));
    D_SG_0   = double(vpa(subs(D_0)));
end