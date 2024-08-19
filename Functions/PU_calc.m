function [Ub,Ib,Zb,Lb,Cb,Ub_dc,Ib_dc,Zb_dc,Cb_dc] = PU_calc(wb,Ul,Sb)

%Base AC-Side:
Ub           = sqrt(2/3)*Ul;
Il           = Sb/(sqrt(3)*Ul);
Ib           = sqrt(2)*Il;
Zb           = Ub / Ib;
Lb           = Zb / wb;
Cb           = 1 / (Zb * wb); 

%Base DC-Side:
Ub_dc        = 2 * Ub;
Ib_dc        = (3/4) * Ib;
Zb_dc        = (8/3) * Zb;
Cb_dc        = (3/8) * Cb;

end