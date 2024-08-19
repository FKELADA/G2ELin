function [A_GFM_0,B_GFM_0,C_GFM_0,D_GFM_0] = GFM_subs(Y_DER,bus_sol,GFM_Data,i_GFM,i_DG,A_0,B_0,C_0,D_0)
%M_AD_0,qm_0,Mvdc_0,pm_0,ved_0,veq_0,igd_0,igq_0,vgd_0,vgq_0,isd_0,isq_0,vmd_0,vmq_0,md_0,mq_0,idc_ref_0,idc_0,M_VLd_0,M_VLq_0,M_CLd_0,M_CLq_0
% GFM Data:
%1)Rf(pu)   2)Lf(pu)   3)Cf(pu)   4)Rt(pu)   5)Lt(pu)  
%6)wb(rad/s)  7)Cdc(pu)  8)Gdc(pu) 9)Kpdc   10)Tdc(s)  
%11)wc(rad/s)  12)wf(rad/s)  13)mp  14)nq  15)Kffi  
%16)Kffv  17)KiVL  18)KpVL  
%19)KiCL 20)KpCL  21)p_ref_0  22)q_ref_0  
%23)ve_ref_0 24)theta_0 25)Eg0_d 26)Eg0_q 27)Ig0_d 28)Ig0_q
%29)Is0_d 30)Is0_q 31)Vm0_d 32)Vm0_q 33)theta_g_0 
%34)dVOC.eta 35)dVOC.alpha 36)VSM.D 37)VSM.J
%38)VSM.Dq 39)VSM.K 40)Matching.K_theta

    Rf        = GFM_Data(i_GFM,1);
    Lf        = GFM_Data(i_GFM,2);
    Cf        = GFM_Data(i_GFM,3);
    Rt        = GFM_Data(i_GFM,4);
    Lt        = GFM_Data(i_GFM,5);
    wff       = 1;
    wb        = GFM_Data(i_GFM,6);
    Cdc       = GFM_Data(i_GFM,7);
    Gdc       = GFM_Data(i_GFM,8);
    Kpdc      = GFM_Data(i_GFM,9);
    Tdc       = GFM_Data(i_GFM,10);
    wc        = GFM_Data(i_GFM,11);
    wf        = GFM_Data(i_GFM,12);
    mp        = GFM_Data(i_GFM,13);
    nq        = GFM_Data(i_GFM,14);
    KiVL      = GFM_Data(i_GFM,17);
    KpVL      = GFM_Data(i_GFM,18);
    Kffi      = GFM_Data(i_GFM,15);
    KiCL      = GFM_Data(i_GFM,19);
    KpCL      = GFM_Data(i_GFM,20);
    Kffv      = GFM_Data(i_GFM,16);
    eta       = GFM_Data(i_GFM,34);
    alfa      = GFM_Data(i_GFM,35);
    Dp        = GFM_Data(i_GFM,36);
    J         = GFM_Data(i_GFM,37);
    Dq        = GFM_Data(i_GFM,38);
    K         = GFM_Data(i_GFM,39);
    K_theta   = GFM_Data(i_GFM,40);

    % Initial Conditions:
    % Setpoints
    p_ref_0   = GFM_Data(i_GFM,21);
    q_ref_0   = GFM_Data(i_GFM,22);
    w_ref_0   = 1;
    ve_ref_0  = GFM_Data(i_GFM,23);
    ved_ref_0 = GFM_Data(i_GFM,23);
    Phi_0     = 0;
    vdc_ref_0 = 1;
    vdc_m_0   = vdc_ref_0;
    
    % Inputs and state variables:
    % States Independent to others:
    M_AD_0    = 0;
    Mvdc_0    = 0;
    qm_0      = q_ref_0;
    pm_0      = p_ref_0;
    vdc_0     = vdc_ref_0;
    w_0       = 1;
    ved_0     = GFM_Data(i_GFM,25);
    veq_0     = GFM_Data(i_GFM,26);
    theta_g_0 = GFM_Data(i_GFM,33);
    % States dependent to others:
    Dw_0      = w_0 - w_ref_0;
    theta_0   = GFM_Data(i_GFM,24);
    vgd_0     = real((bus_sol(Y_DER(i_DG,1),2))*exp(1i*((bus_sol(Y_DER(i_DG,1),3))*pi/180-theta_0)));
    vgq_0     = imag((bus_sol(Y_DER(i_DG,1),2))*exp(1i*((bus_sol(Y_DER(i_DG,1),3))*pi/180-theta_0)));
    vgd_g_0   = vgd_0*cos(theta_g_0 - theta_0) + vgq_0*sin(theta_g_0 - theta_0);
    vgq_g_0   = -vgd_0*sin(theta_g_0 - theta_0) + vgq_0*cos(theta_g_0 - theta_0);
    igd_0     = GFM_Data(i_GFM,27); 
    igq_0     = GFM_Data(i_GFM,28); 
    isd_0     = GFM_Data(i_GFM,29);
    isq_0     = GFM_Data(i_GFM,30);
    vmd_0     = GFM_Data(i_GFM,31);
    vmq_0     = GFM_Data(i_GFM,32);
    md_0      = vmd_0/vdc_ref_0;
    mq_0      = vmq_0/vdc_ref_0;
    idc_ref_0 = p_ref_0/vdc_ref_0;
    idc_0     = idc_ref_0;
    M_VLd_0   = (isd_0 - igd_0*Kffi + veq_0*Cf*wff);
    M_VLq_0   = (isq_0 - igq_0*Kffi - ved_0*Cf*wff);
    M_CLd_0   = (vmd_0 - ved_0*Kffv + isq_0*Lf*wff);
    M_CLq_0   = (vmq_0 - veq_0*Kffv - isd_0*Lf*wff);

    A_GFM_0   = double(vpa(subs(A_0)));
    B_GFM_0   = double(vpa(subs(B_0)));
    C_GFM_0   = double(vpa(subs(C_0)));
    D_GFM_0   = double(vpa(subs(D_0)));
end