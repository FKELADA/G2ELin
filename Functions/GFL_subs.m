function [A_GFL_0,B_GFL_0,C_GFL_0,D_GFL_0] = GFL_subs(Y_DER,bus_sol,GFL_Data,i_GFL,i_DG,A_0,B_0,C_0,D_0)
% GFL Data:
%1)Rf(pu)   2)Lf(pu)   3)Cf(pu)   4)Rt(pu)   5)Lt(pu)  
%6)wb(rad/s)  7)Cdc(pu)  8)Gdc(pu)   9)Tdc(s)  
%10)Kpd  11)Kid  12)Kiq  13)Kffv  14)KpPLL  15)KiPLL  
%16)KiCL 17)KpCL  18)p_ref_0  19)q_ref_0  
%20)theta_0 21)Eg0_d 22)Eg0_q 23)Ig0_d 24)Ig0_q
%25)Is0_d 26)Is0_q 27)Vm0_d 28)Vm0_q 29)M_pll_0
%30)M_d_0 31)M_q_0 32)M_CLd_0 33)M_CLq_0

    Rf        = GFL_Data(i_GFL,1);
    Lf        = GFL_Data(i_GFL,2);
    Cf        = GFL_Data(i_GFL,3);
    Rt        = GFL_Data(i_GFL,4);
    Lt        = GFL_Data(i_GFL,5);
    wff       = 1;
    wb        = GFL_Data(i_GFL,6);
    Cdc       = GFL_Data(i_GFL,7);
    Gdc       = GFL_Data(i_GFL,8);
    Tdc       = GFL_Data(i_GFL,9);
    Kpd       = GFL_Data(i_GFL,10);
    Kid       = GFL_Data(i_GFL,11);
    Kiq       = GFL_Data(i_GFL,12);
    Kffv      = GFL_Data(i_GFL,13);
    Kppll     = GFL_Data(i_GFL,14);
    Kipll     = GFL_Data(i_GFL,15);
    KiCL      = GFL_Data(i_GFL,16);
    KpCL      = GFL_Data(i_GFL,17);

    % Initial Conditions:
    % Setpoints
    p_ref_0   = GFL_Data(i_GFL,18);
    q_ref_0   = GFL_Data(i_GFL,19);
    vdc_ref_0 = 1;
    idc_ref_0 = p_ref_0/vdc_ref_0;
    vdc_0     = vdc_ref_0;
    theta_g_0   = GFL_Data(i_GFL,34);
    theta_pll_0 = GFL_Data(i_GFL,20);
    w_pll_0   = 1;
    vgd_0     = real((bus_sol(Y_DER(i_DG,1),2))*exp(1i*((bus_sol(Y_DER(i_DG,1),3))*pi/180-theta_pll_0)));
    vgq_0     = imag((bus_sol(Y_DER(i_DG,1),2))*exp(1i*((bus_sol(Y_DER(i_DG,1),3))*pi/180-theta_pll_0)));
    vgd_g_0   = vgd_0*cos(theta_g_0 - theta_pll_0) + vgq_0*sin(theta_g_0 - theta_pll_0);
    vgq_g_0   = -vgd_0*sin(theta_g_0 - theta_pll_0) + vgq_0*cos(theta_g_0 - theta_pll_0);
    
    % Inputs and state variables:
    % States Independent to others:
    qm_0      = q_ref_0;
    ved_0     = GFL_Data(i_GFL,21);
    veq_0     = GFL_Data(i_GFL,22);

    
    % States dependent to others:
    igd_0     = GFL_Data(i_GFL,23); 
    igq_0     = GFL_Data(i_GFL,24); 
    M_pll_0   = GFL_Data(i_GFL,29);
    isd_0     = GFL_Data(i_GFL,25);
    M_d_0     = GFL_Data(i_GFL,30);
    isq_0     = GFL_Data(i_GFL,26);
    M_q_0     = GFL_Data(i_GFL,31);
    vmd_0     = GFL_Data(i_GFL,27);
    vmq_0     = GFL_Data(i_GFL,28);
    md_0      = vmd_0/vdc_ref_0;
    mq_0      = vmq_0/vdc_ref_0;
    idc_0     = idc_ref_0;
    M_CLd_0   = GFL_Data(i_GFL,32);
    M_CLq_0   = GFL_Data(i_GFL,33);

    A_GFL_0   = double(vpa(subs(A_0)));
    B_GFL_0   = double(vpa(subs(B_0)));
    C_GFL_0   = double(vpa(subs(C_0)));
    D_GFL_0   = double(vpa(subs(D_0)));
end