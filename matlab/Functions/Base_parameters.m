function [IB_BP,SM_BP,TR1_BP,TR2_BP,Lines_BP,Loads_BP,VSC_BP] = Base_parameters(Voltage_level,wb,Sb,Un_LV,Un_MV,Un_HV)

[Ub_LV,Ib_LV,Zb_LV,Lb_LV,Cb_LV,Ub_LV_dc,Ib_LV_dc,Zb_LV_dc,Cb_LV_dc] = PU_calc(wb,Un_LV,Sb);
[Ub_MV,Ib_MV,Zb_MV,Lb_MV,Cb_MV,Ub_MV_dc,Ib_MV_dc,Zb_MV_dc,Cb_MV_dc] = PU_calc(wb,Un_MV,Sb);
[Ub_HV,Ib_HV,Zb_HV,Lb_HV,Cb_HV,Ub_HV_dc,Ib_HV_dc,Zb_HV_dc,Cb_HV_dc] = PU_calc(wb,Un_HV,Sb);

switch Voltage_level(1)
    case 'HV'
        IB.Un     = Un_HV;
        IB.Ub     = Ub_HV;
        IB.Ib     = Ib_HV;
        IB.Zb     = Zb_HV;
        IB.Lb     = Lb_HV;
        IB.Cb     = Cb_HV;
    case 'MV'
        IB.Un     = Un_MV;
        IB.Ub     = Ub_MV;
        IB.Ib     = Ib_MV;
        IB.Zb     = Zb_MV;
        IB.Lb     = Lb_MV;
        IB.Cb     = Cb_MV;        
    case 'LV'
        IB.Un     = Un_LV;
        IB.Ub     = Ub_LV;
        IB.Ib     = Ib_LV;
        IB.Zb     = Zb_LV;
        IB.Lb     = Lb_LV;
        IB.Cb     = Cb_LV;        
end

IB_BP = [IB.Un IB.Ub IB.Ib IB.Zb IB.Lb IB.Cb];

switch Voltage_level(2)
    case 'HV'
        SM.Un     = Un_HV;
        SM.Es_base= Ub_HV;
        SM.Is_base= Ib_HV;
        SM.Cb     = Cb_HV;
        SM.Zb     = Zb_HV;
    case 'MV'
        SM.Un     = Un_MV;
        SM.Es_base= Ub_MV;
        SM.Is_base= Ib_MV;
        SM.Cb     = Cb_MV;
        SM.Zb     = Zb_MV;
    case 'LV'
        SM.Un     = Un_LV;
        SM.Es_base= Ub_LV;
        SM.Is_base= Ib_LV;
        SM.Cb     = Cb_LV;
        SM.Zb     = Zb_LV;
end

SM_BP = [SM.Un SM.Es_base SM.Is_base SM.Cb SM.Zb];

% for the TRs primary:
switch Voltage_level(3)
    case 'HV'
        TR1.Un     = Un_HV;
        TR1.Ub     = Ub_HV;
        TR1.Ib     = Ib_HV;
        TR1.Zb     = Zb_HV;
        TR1.Lb     = Lb_HV;
    case 'MV'
        TR1.Un     = Un_MV;
        TR1.Ub     = Ub_MV;
        TR1.Ib     = Ib_MV;
        TR1.Zb     = Zb_MV;
        TR1.Lb     = Lb_MV;      
    case 'LV'
        TR1.Un     = Un_LV;
        TR1.Ub     = Ub_LV;
        TR1.Ib     = Ib_LV;
        TR1.Zb     = Zb_LV;
        TR1.Lb     = Lb_LV;     
end

TR1_BP = [TR1.Un TR1.Ub TR1.Ib TR1.Zb TR1.Lb];

% for the TRs seondary:
switch Voltage_level(4)
    case 'HV'
        TR2.Un     = Un_HV;
        TR2.Ub     = Ub_HV;
        TR2.Ib     = Ib_HV;
        TR2.Zb     = Zb_HV;
        TR2.Lb     = Lb_HV;
    case 'MV'
        TR2.Un     = Un_MV;
        TR2.Ub     = Ub_MV;
        TR2.Ib     = Ib_MV;
        TR2.Zb     = Zb_MV;
        TR2.Lb     = Lb_MV;      
    case 'LV'
        TR2.Un     = Un_LV;
        TR2.Ub     = Ub_LV;
        TR2.Ib     = Ib_LV;
        TR2.Zb     = Zb_LV;
        TR2.Lb     = Lb_LV;     
end

TR2_BP = [TR2.Un TR2.Ub TR2.Ib TR2.Zb TR2.Lb];

% for the Lines:
switch Voltage_level(5)
    case 'HV'
        Line.Un     = Un_HV;
        Line.Ub     = Ub_HV;
        Line.Ib     = Ib_HV;
        Line.Zb     = Zb_HV;
        Line.Lb     = Lb_HV;
        Line.Cb     = Cb_HV;
    case 'MV'
        Line.Un     = Un_MV;
        Line.Ub     = Ub_MV;
        Line.Ib     = Ib_MV;
        Line.Zb     = Zb_MV;
        Line.Lb     = Lb_MV;
        Line.Cb     = Cb_MV;        
    case 'LV'
        Line.Un     = Un_LV;
        Line.Ub     = Ub_LV;
        Line.Ib     = Ib_LV;
        Line.Zb     = Zb_LV;
        Line.Lb     = Lb_LV;
        Line.Cb     = Cb_LV;        
end

Lines_BP = [Line.Un Line.Ub Line.Ib Line.Zb Line.Lb Line.Cb];

% for the Loads:
switch Voltage_level(6)
    case 'HV'
        Load.Un     = Un_HV;
        Load.Ub     = Ub_HV;
        Load.Ib     = Ib_HV;
        Load.Zb     = Zb_HV;
        Load.Lb     = Lb_HV;
        Load.Cb     = Cb_HV;
    case 'MV'
        Load.Un     = Un_MV;
        Load.Ub     = Ub_MV;
        Load.Ib     = Ib_MV;
        Load.Zb     = Zb_MV;
        Load.Lb     = Lb_MV;
        Load.Cb     = Cb_MV;        
    case 'LV'
        Load.Un     = Un_LV;
        Load.Ub     = Ub_LV;
        Load.Ib     = Ib_LV;
        Load.Zb     = Zb_LV;
        Load.Lb     = Lb_LV;
        Load.Cb     = Cb_LV;        
end

Loads_BP = [Load.Un Load.Ub Load.Ib Load.Zb Load.Lb Load.Cb];

% for the VSCS:
switch Voltage_level(7)
    case 'HV'
        VSC.Un     = Un_HV;
        VSC.Ub     = Ub_HV;
        VSC.Ib     = Ib_HV;
        VSC.Zb     = Zb_HV;
        VSC.Lb     = Lb_HV;
        VSC.Cb     = Cb_HV;
        VSC.Zb_dc  = Zb_HV_dc;
        VSC.Ub_dc  = Ub_HV_dc;
        VSC.Ib_dc  = Ib_HV_dc;
        VSC.Cb_dc  = Cb_HV_dc;
    case 'MV'
        VSC.Un     = Un_MV;
        VSC.Ub     = Ub_MV;
        VSC.Ib     = Ib_MV;
        VSC.Zb     = Zb_MV;
        VSC.Lb     = Lb_MV;
        VSC.Cb     = Cb_MV;        
        VSC.Ub_dc  = Ub_MV_dc;
        VSC.Ib_dc  = Ib_MV_dc;
        VSC.Zb_dc  = Zb_MV_dc;
        VSC.Cb_dc  = Cb_MV_dc;
    case 'LV'
        VSC.Un     = Un_LV;
        VSC.Ub     = Ub_LV;
        VSC.Ib     = Ib_LV;
        VSC.Zb     = Zb_LV;
        VSC.Lb     = Lb_LV;
        VSC.Cb     = Cb_LV;        
        VSC.Ub_dc  = Ub_LV_dc;
        VSC.Ib_dc  = Ib_LV_dc;
        VSC.Zb_dc  = Zb_LV_dc;
        VSC.Cb_dc  = Cb_LV_dc;
end

VSC_BP = [VSC.Un VSC.Ub VSC.Ib VSC.Zb VSC.Lb VSC.Cb VSC.Ub_dc VSC.Ib_dc VSC.Zb_dc VSC.Cb_dc];

end