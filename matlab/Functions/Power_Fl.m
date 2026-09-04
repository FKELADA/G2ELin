function [bus_sol,P_loss,I_lignes_pu,iter,conv_flag,Y] = Power_Fl(line,bus,Sb_MW) 

% Inputs parameters : matrices that describe the network 

% Line : based on the Pi model
%            d        R       X=lw      a
%       -------------[]------mmm-------------
%            |                          |
%          -----                      -----  Cw/2
%          -----                      -----
%            |                          |
% 
% - Column 1 : Departure node
% - Column 2 : Arriving node
% - Column 3 : R (pu)
% - Column 4 : X(pu)
% - Column 5 : Cw/2 (pu)


% Bus
%   - Column 1 : Number of node
%   - Column 2 : Initial voltage at the nodes
%   - Column 3 : Initial voltage phase angle regarding reference node (slack bus)
%   - Column 4 : Active power generated (en pu)
%   - Column 5 : Reactive power generated (en pu)
%   - Column 6 : Active power consummed (en pu)
%   - Column 7 : Reactive power consummed (en pu)
%   - Column 8 : not used here (zeros)
%   - Column 9 : not used here (zeros)
%   - Column 10 : Type of node (1 for slack, 2 for PV and 3 for PQ)

Vmin = 0.95;
Vmax = 1.05;

% Outputs parameters
% - bus_sol : solution bus matrix : 
%   - Column 1 : Number of node
%   - Column 2 : Final voltage at the nodes computed by the loadflow
%   - Column 3 : Final voltage phase angle regarding reference node (slack bus) computed by the loadflow
%   (noeud bilan)au noeud calculée par le loadflow
%   - Column 4 : Active power generated (en pu)
%   - Column 5 : Rective power generated (en pu)
%   - Column 6 : Active power consummed (en pu)
%   - Column 7 : Reactive power consummed (en pu)
%   - Column 8 : not used here (zeros)
%   - Column 9 : not used here (zeros)
%   - Column 10 : Type of node
%   - Column 10 : not used here 
%   - Column 11 : not used here
%   - Column 12 : not used here (ones)
%   - Column 13 : Vmax 
%   - Column 14 : Vmin 

% - P_loss : power losses computed by the loadflow (in pu)

% - I_lignes_pu : Current in lines (in pu)

% - iter = number of iteration of the loadflow




    global bus_int
    global Qg bus_type g_bno PQV_no PQ_no ang_red volt_red 
    global Q Ql
    global gen_chg_idx
    nbus = length(bus(:,1));
    gen_chg_idx = ones(nbus,1);
    tol=1e-9;
    iter_max=20;
    acc=1.0;
    display='y';
    flag=2;
    %tap=1;

    tt = clock;     % start the total time clock
    jay = sqrt(-1);
    load_bus = 3;
    gen_bus = 2;
    swing_bus = 1;
    if exist('flag') == 0
      flag = 1;
    end
    lf_flag = 1;
    % set solution defaults
    if isempty(tol); tol = 1e-9;end
    if isempty(iter_max);iter_max = 30;end
    if isempty(acc);acc = 1.0; end;
    if isempty(display);display = 'n';end;

    if flag <1 | flag > 2
      error('LOADFLOW: flag not recognized')
    end
    [nline nlc] = size(line);    % number of lines and no of line cols
    [nbus ncol] = size(bus);     % number of buses and number of col
    % set defaults
    % bus data defaults
    if ncol<15
        % set generator var limits
        if ncol<12
          bus(:,11) = 9999*ones(nbus,1);
          bus(:,12) = -9999*ones(nbus,1);
        end 
        if ncol<13;bus(:,13) = ones(nbus,1);end
        bus(:,14) = Vmax*ones(nbus,1);
        bus(:,15) = Vmin*ones(nbus,1);
        volt_min = bus(:,15);
        volt_max = bus(:,14);
    else
       volt_min = bus(:,15);
       volt_max = bus(:,14);
    end
    no_vmin_idx = find(volt_min==0);
    if ~isempty(no_vmin_idx)
      volt_min(no_vmin_idx) = Vmin*ones(length(no_vmin_idx),1);
    end
    no_vmax_idx = find(volt_max==0);
    if ~isempty(no_vmax_idx)
      volt_max(no_vmax_idx) = Vmax*ones(length(no_vmax_idx),1);
    end
    no_mxv = find(bus(:,11)==0);
    no_mnv = find(bus(:,12)==0);
    if ~isempty(no_mxv);bus(no_mxv,11)=9999*ones(length(no_mxv),1);end
    if ~isempty(no_mnv);bus(no_mnv,12) = -9999*ones(length(no_mnv),1);end
    no_vrate = find(bus(:,13)==0);
    if ~isempty(no_vrate);bus(no_vrate,13) = ones(length(no_vrate),1);end
    tap_it = 0;
    tap_it_max = 30;
    no_taps = 0;
    % line data defaults, sets all tap ranges to zero - this fixes taps
    if nlc < 10
      line(:,7:10) = zeros(nline,4);
      no_taps = 1;
      % disable tap changing
    end

    % outer loop for on-load tap changers

    mm_chk=1;
    while (tap_it<tap_it_max&mm_chk)
      tap_it = tap_it+1;
      % build admittance matrix Y
      [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line);
      % process bus data
      bus_no = bus(:,1);
      V = bus(:,2);
      ang = bus(:,3)*pi/180;
      Pg = bus(:,4);
      Qg = bus(:,5);
      Pl = bus(:,6);
      Ql = bus(:,7);
      Gb = bus(:,8);
      Bb = bus(:,9);
      bus_type = round(bus(:,10));
      qg_max = bus(:,11);
      qg_min = bus(:,12);
      sw_bno=ones(nbus,1);
      g_bno=sw_bno;
      % set up index for Jacobian calculation
      %% form PQV_no and PQ_no
      bus_zeros=zeros(nbus,1);
      swing_index=find(bus_type==1); 
      sw_bno(swing_index)=bus_zeros(swing_index);
      PQV_no=find(bus_type >=2);
      PQ_no=find(bus_type==3);
      gen_index=find(bus_type==2);
      g_bno(gen_index)=bus_zeros(gen_index); 
      %sw_bno is a vector having ones everywhere but the swing bus locations
      %g_bno is a vector having ones everywhere but the generator bus locations

      % construct sparse angle reduction matrix
      il = length(PQV_no);
      ii = [1:1:il]';
      ang_red = sparse(ii,PQV_no,ones(il,1),il,nbus);

      % construct sparse voltage reduction matrix
      il = length(PQ_no);
      ii = [1:1:il]';
      volt_red = sparse(ii,PQ_no,ones(il,1),il,nbus);

      iter = 0;     % initialize iteration counter

    % calculate the power mismatch and check convergence

      [delP,delQ,P,Q,conv_flag,mism] =...
                 calc(nbus,V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol);


      st = clock;     % start the iteration time clock
      %% start iteration process for main Newton_Raphson solution
      while (conv_flag == 1 & iter < iter_max)
         iter = iter + 1;
         % Form the Jacobean matrix
         clear Jac
         Jac=form_jac(V,ang,Y,ang_red,volt_red);
         % reduced real and reactive power mismatch vectors
         red_delP = ang_red*delP;
         red_delQ = volt_red*delQ;
         clear delP delQ
         % solve for voltage magnitude and phase angle increments
         temp = Jac\[red_delP; red_delQ];
         % expand solution vectors to all buses
         delAng = ang_red'*temp(1:length(PQV_no),:);
         delV = volt_red'*temp(length(PQV_no)+1:length(PQV_no)+length(PQ_no),:);
         % update voltage magnitude and phase angle
         V = V + acc*delV;
         %V = max(V,volt_min);  % voltage higher than minimum
         %V = min(V,volt_max);  % voltage lower than maximum
         ang = ang + acc*delAng;
         % calculate the power mismatch and check convergence
         [delP,delQ,P,Q,conv_flag,mism] =...
                   calc(nbus,V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol);
         % check if Qg is outside limits
         gen_index=find(bus_type==2);
         Qg(gen_index) = Q(gen_index) + Ql(gen_index);
         lim_flag = chq_lim(qg_max,qg_min);
         if lim_flag == 1; 
           disp('Qg at var limit');
         end
      end
      if no_taps == 0
    %%%%%%%
    %lftap%
    %%%%%%%
       max_v_idx = find(V>=volt_max);
       min_v_idx = find(V<=volt_min);
       mm_chk = 0;
       if (~isempty(max_v_idx)||~isempty(min_v_idx))
         mm_chk = 1;
         if ~isempty(max_v_idx)
         % change on load taps to correct voltage
         % assumes that bus to be corrected is the to bus
           for fb = 1 : length(max_v_idx)
             chk_fb = find(line(:,2) == bus(max_v_idx(fb),1));
             if any(line(chk_fb,8)~=0),
               if line(chk_fb(find(line(chk_fb,8)~=0)),8)~=0
                 % can change tap
                 % get from bus voltage
    %              disp('voltage high changing tap on line');disp(chk_fb(find(line(chk_fb,8)~=0)))
                 vm1 = V(bus_int(line(chk_fb(find(line(chk_fb,8)~=0)),1)));
                 vm2 = volt_max(max_v_idx(fb));
                 verror = V(max_v_idx(fb))-vm2;
                 % voltage too high, tap must be increased
                 if verror/vm2<line(chk_fb(find(line(chk_fb,8)~=0)),10)
                    % increase tap by one step
                    tap = line(chk_fb(find(line(chk_fb,8)~=0)),6) + line(chk_fb(find(line(chk_fb,8)~=0)),10);
                 else
                    tap = line(chk_fb(find(line(chk_fb,8)~=0)),6) + verror/vm2;
                    tap_set = ceil( (tap-line(chk_fb(find(line(chk_fb,8)~=0)),9))/line(chk_fb(find(line(chk_fb,8)~=0)),10));
                    tap = tap_set*line(chk_fb(find(line(chk_fb,8)~=0)),10) + line(chk_fb(find(line(chk_fb,8)~=0)),9);
                 end
                 tap = min(line(chk_fb(find(line(chk_fb,8)~=0)),8),max(tap, line(chk_fb(find(line(chk_fb,8)~=0)),9)));
                 % reset tap in line data
    %              disp('tap reset to');tap
                 line(chk_fb(find(line(chk_fb,8)~=0)),6) = tap; %#ok<*FNDSB>
               end
             end
           end
         end 
         if ~isempty(min_v_idx)
         % change on load taps to correct voltage
         % assumes that bus to be corrected is the to bus
           for fb = 1 : length(min_v_idx)
             chk_fb = find(line(:,2) == bus(min_v_idx(fb),1));
             if ~isempty(chk_fb)
               if any(line(chk_fb,8)~=0),
                 % can change tap
    %              disp('voltage low changing tap on line');disp(chk_fb(find(line(chk_fb,8)~=0)))
                 vm1 = V(bus_int(line(chk_fb(find(line(chk_fb,8)~=0)),1)));
                 vm2 = volt_min(min_v_idx(fb));
                 verror = vm2 - V(min_v_idx(fb));
                 % voltage too low tap must be reduced
                 if verror/vm2<line(chk_fb(find(line(chk_fb,8)~=0)),10)
                   % reduce tap by one increment
                   tap = line(chk_fb(find(line(chk_fb,8)~=0)),6)-line(chk_fb(find(line(chk_fb,8)~=0)),10);
                 else
                   tap = line(chk_fb(find(line(chk_fb,8)~=0)),6) - verror./vm2;
                   tap_set = fix( (tap-line(chk_fb(find(line(chk_fb,8)~=0)),9))/line(chk_fb(find(line(chk_fb,8)~=0)),10));
                   tap = tap_set*line(chk_fb(find(line(chk_fb,8)~=0)),10) + line(chk_fb(find(line(chk_fb,8)~=0)),9);
                 end
                 tap = min(line(chk_fb(find(line(chk_fb,8)~=0)),8),max(tap, line(chk_fb(find(line(chk_fb,8)~=0)),9)));
    %              disp('taps reset to');tap
                 % reset tap in line data
                 line(chk_fb(find(line(chk_fb,8)~=0)),6) = tap;
               end
             end
           end
         end
       end
       %%%%%%%%
      else
         mm_chk = 0;
      end
    end
    % if tap_it >= tap_it_max
    % %   disp('tap iteration failed to converge after 10 iterations')
    % end
    ste = clock;     % end the iteration time clock
    vmx_idx = find(V==volt_max);
    vmn_idx = find(V==volt_min);
    if ~isempty(vmx_idx)
%       disp('voltages at') 
%       bus(vmx_idx,1)' 
%       disp('are at the max limit')
    end
    if ~isempty(vmn_idx)
%       disp('voltages at') 
%       bus(vmn_idx,1)' 
%       disp('are at the min limit');
    end 
    gen_index=find(bus_type==2);
    load_index = find(bus_type==3);
    Pg(gen_index) = P(gen_index) + Pl(gen_index);
    Qg(gen_index) = Q(gen_index) + Ql(gen_index);
    gend_idx = find((bus(:,10)==2)&(bus_type~=2));
    if ~isempty(gend_idx)
       Qg(gend_idx) = bus(gend_idx,7)-Ql(gend_idx);
       disp('the following generators are at their var limits')
       disp('    bus#    Qg')
       disp([bus(gend_idx,1)  Qg(gend_idx)])
       Ql(gend_idx)=Ql(gend_idx)-Qg(gend_idx);
       Qg(gend_idx)=zeros(length(gend_idx),1);
    end
    Pl(load_index) = Pg(load_index) - P(load_index);
    Ql(load_index) = Qg(load_index) - Q(load_index);

    Pg(SB) = P(SB) + Pl(SB); Qg(SB) = Q(SB) + Ql(SB);
    VV = V.*exp(jay*ang);  % solution voltage 
    % calculate the line flows and power losses
    tap_index = find(abs(line(:,6))>0);
    tap_ratio = ones(nline,1);
    tap_ratio(tap_index)=line(tap_index,6);
    phase_shift(:,1) = line(:,7);
    tps = tap_ratio.*exp(jay*phase_shift*pi/180);
    from_bus = line(:,1);
    from_int = bus_int(round(from_bus));
    to_bus = line(:,2);
    to_int = bus_int(round(to_bus));
    r = line(:,3);
    rx = line(:,4);
    chrg = line(:,5);
    z = r + jay*rx;
    y = ones(nline,1)./z;


    MW_s = VV(from_int).*conj((VV(from_int) - tps.*VV(to_int)).*y ...
           + VV(from_int).*(jay*chrg/2))./(tps.*conj(tps));
    P_s = real(MW_s);     % active power sent out by from_bus
                          % to to_bus
    Q_s = imag(MW_s);     % reactive power sent out by 
                          % from_bus to to_bus
    MW_r = VV(to_int).*conj((VV(to_int) ...
           - VV(from_int)./tps).*y ...
           + VV(to_int).*(jay*chrg/2));
    P_r = real(MW_r);     % active power received by to_bus 
                          % from from_bus
    Q_r = imag(MW_r);     % reactive power received by 
                          % to_bus from from_bus
    iline = [1:1:nline]';
      line_ffrom = [iline from_bus to_bus P_s Q_s];
      line_fto   = [iline to_bus from_bus P_r Q_r];
    % keyboard
    P_loss = sum(P_s) + sum(P_r) ;
    Q_loss = sum(Q_s) + sum(Q_r) ;
    bus_sol=[bus_no  V  ang*180/pi Pg Qg Pl Ql Gb Bb...
             bus_type qg_max qg_min bus(:,13) volt_max volt_min];
    line_sol = line;
    line_flow(1:nline, :)  =[iline from_bus to_bus P_s Q_s];
    line_flow(1+nline:2*nline,:) = [iline to_bus from_bus P_r Q_r]; 
    % display results
    if display == 'y',
     clc
     disp('                             LOAD-FLOW STUDY')
     disp('                    REPORT OF POWER FLOW CALCULATIONS ')
     disp(' ')
     disp(date)
     fprintf('SWING BUS                  : BUS %g \n', SB)
        fprintf('NUMBER OF ITERATIONS       : %g \n', iter)
        fprintf('SOLUTION TIME              : %g sec.\n',etime(ste,st))
        fprintf('TOTAL TIME                 : %g sec.\n',etime(clock,tt))
        fprintf('TOTAL REAL POWER GENERATION    : %g MW.\n',sum(bus_sol(:,4))*Sb_MW)
        fprintf('TOTAL REACTIVE POWER GENERATION    : %g MVAr.\n',sum(bus_sol(:,5))*Sb_MW)
        fprintf('TOTAL REAL POWER LOAD    : %g MW.\n',sum(bus_sol(:,6))*Sb_MW)
        fprintf('TOTAL REAL REACTIVE LOAD    : %g MVAr.\n',sum(bus_sol(:,7))*Sb_MW)  
        fprintf('TOTAL REAL POWER LOSSES    : %g MW.\n',P_loss*Sb_MW)
        fprintf('TOTAL REACTIVE POWER LOSSES: %g MVAr.\n\n',Q_loss*Sb_MW)
      if conv_flag == 0,
       disp('                                      GENERATION             LOAD')
       disp('       BUS     VOLTS     ANGLE(°)   REAL  REACTIVE      REAL  REACTIVE ')
       disp(bus_sol(:,1:7))

       disp('                      LINE FLOWS                     ')
       disp('      LINE  FROM BUS    TO BUS      REAL  REACTIVE   ')
       disp(line_ffrom)
       disp(line_fto)
      end
    end; %
    if iter > iter_max,
      disp('Note: Solution did not converge in %g iterations.\n', iter_max)
      lf_flag = 0
    end

I_lignes_pu=abs(((VV(from_int) - tps.*VV(to_int)).*y ...
           + VV(from_int).*(jay*chrg/2))./(tps.*conj(tps)));

function [delP,delQ,P,Q,conv_flag,mism] = ...
                 calc(nbus,V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol)
    % Syntax:  [delP,delQ,P,Q,conv_flag] = 
    %                calc(nbus,V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol)
    %
    % Purpose: calculates power mismatch and checks convergence
    %          also determines the values of P and Q based on the 
    %          supplied values of voltage magnitude and angle
    % Version: 2.0 eliminates do loop
    % Input:   nbus      - total number of buses
    %          bus_type  - load_bus(3), gen_bus(2), swing_bus(1)
    %          V         - magnitude of bus voltage
    %          ang       - angle(rad) of bus voltage
    %          Y         - admittance matrix
    %          Pg        - real power of generation
    %          Qg        - reactive power  of generation
    %          Pl        - real power of load
    %          Ql        - reactive power of load
    %	 sw_bno - a vector having zeros at all  swing_bus locations ones otherwise
    %	 g_bno  - a vector having zeros at all  generator bus locations ones otherwise
    %          tol       - a tolerance of computational error
    %
    % Output:  delP      - real power mismatch
    %          delQ      - reactive power mismatch
    %          P         - calculated real power
    %          Q         - calculated reactive power
    %          conv_flag - 0, converged
    %                      1, not yet converged
    %
    % See also:  
    %
    % Calls:
    %
    % Called By:   loadflow

    % (c) Copyright 1991 Joe H. Chow - All Rights Reserved
    %
    % History (in reverse chronological order)
    % Version:   2.0
    % Author:    Graham Rogers
    % Date:      July 1994
    %
    % Version:   1.0
    % Author:    Kwok W. Cheung, Joe H. Chow
    % Date:      March 1991
    %
    % ************************************************************
    jay = sqrt(-1);
    swing_bus = 1;
    gen_bus = 2;
    load_bus = 3;
    % voltage in rectangular coordinate
    V_rect = V.*exp(jay*ang);  
    % bus current injection
    cur_inj = Y*V_rect;
    % power output based on voltages 
    S = V_rect.*conj(cur_inj);
    P = real(S); Q = imag(S);
    delP = Pg - Pl - P;
    delQ = Qg - Ql - Q;
    % zero out mismatches on swing bus and generation bus
    delP=delP.*sw_bno;
    delQ=delQ.*sw_bno;
    delQ=delQ.*g_bno;
    %  total mismatch
    [pmis,ip]=max(abs(delP));
    [qmis,iq]=max(abs(delQ));
    mism = pmis+qmis;
    if mism > tol,
        conv_flag = 1;
      else
        conv_flag = 0;
    end
    return

function f = chq_lim(qg_max,qg_min)
%Syntax:
%       f = chq_lim(qg_max,qg_min)
% function for detecting generator vars outside limit
% sets Qg to zero if limit exceded, sets Ql to negative of limit
% sets bus_type to 3, and recalculates ang_red and volt_red
% changes  generator bus_type to type 3
% recalculates the generator index
% inputs: qg_max and qg_min are the last two clumns of the bus matrix
% outputs:f is set to zero if no limit reached, or to 1 if a limit is reached
% Version:  1.1
% Author: Graham Rogers
% Date:   May 1997
% Purpose: Addition of var limit index
% Version:  1.0
% Author:   Graham Rogers
% Date:     October 1996
%
% (c) copyright Joe Chow 1996
global Qg bus_type g_bno PQV_no PQ_no ang_red volt_red 
global Q Ql
global gen_chg_idx
%         gen_chg_idx indicates those generators changed to PQ buses
%         gen_cgq_idx = ones(n of bus,1) if no gen at vars limits
%                     = 0 at the corresponding bus if generator at var limit 

f = 0;
lim_flag = 0;% indicates whether limit has been reached
gen_idx = find(bus_type ==2);
qg_max_idx = find(Qg(gen_idx)>qg_max(gen_idx));
qg_min_idx = find(Qg(gen_idx)<qg_min(gen_idx));
if ~isempty(qg_max_idx)
  %some q excedes maximum
  %set Qg to zero
  Qg(gen_idx(qg_max_idx)) = zeros(length(qg_max_idx),1);
  % modify Ql
  Ql(gen_idx(qg_max_idx)) = Ql(gen_idx(qg_max_idx))...
                            - qg_max(gen_idx(qg_max_idx));
  % modify bus_type to PQ bus
  bus_type(gen_idx(qg_max_idx)) = 3*ones(length(qg_max_idx),1);
  gen_chg_idx(gen_idx(qg_max_idx)) = zeros(length(qg_max_idx),1);
  lim_flag = 1;
end
if ~isempty(qg_min_idx)
  %some q less than minimum
  %set Qg to zero
  Qg(gen_idx(qg_min_idx)) = zeros(length(qg_min_idx),1);
  % modify Ql
  Ql(gen_idx(qg_min_idx)) = Ql(gen_idx(qg_min_idx))...
                            - qg_min(gen_idx(qg_min_idx));
  % modify bus_type to PQ bus
  bus_type(gen_idx(qg_min_idx)) = 3*ones(length(qg_min_idx),1);
  gen_chg_idx(gen_idx(qg_min_idx)) = zeros(length(qg_min_idx),1);
  lim_flag = 1;
end
if lim_flag == 1
  %recalculate g_bno
  nbus = length(bus_type);
  g_bno = ones(nbus,1);
  bus_zeros=zeros(nbus,1);
  bus_index=[1:1:nbus]';
  PQV_no=find(bus_type >=2);
  PQ_no=find(bus_type==3);
  gen_index=find(bus_type==2);
  g_bno(gen_index)=bus_zeros(gen_index); 
  % construct sparse angle reduction matrix
  il = length(PQV_no);
  ii = [1:1:il]';
  ang_red = sparse(ii,PQV_no,ones(il,1),il,nbus);
  % construct sparse voltage reduction matrix
  il = length(PQ_no);
  ii = [1:1:il]';
  volt_red = sparse(ii,PQ_no,ones(il,1),il,nbus);
end
f = lim_flag;
return

function [Jac11,Jac12,Jac21,Jac22]=form_jac(V,ang,Y,ang_red,volt_red)
    % Syntax:  [Jac] = form_jac(V,ang,Y,ang_red,volt_red)
    %          [Jac11,Jac12,Jac21,Jac22] = form_jac(V,ang,Y,...
    %                                      ang_red,volt_red)
    %
    % Purpose: form the Jacobian matrix using sparse matrix techniques
    %
    % Input:   V        - magnitude of bus voltage
    %          ang      - angle(rad) of bus voltage
    %          Y        - admittance matrix
    %          ang_red  - matrix to eliminate swing bus  voltage magnitude and angle 
    %                     entries
    %          volt_red - matrix to eliminate generator bus voltage magnitude
    %                     entries
    % Output:  Jac      - jacobian matrix
    %          Jac11,Jac12,Jac21,Jac22 - submatrices of 
    %                                      jacobian matrix  
    % See also:   
    %
    % Calls:
    %
    % Called By:   vsdemo loadflow

    % (c) Copyright 1991-1996 Joe H. Chow - All Rights Reserved
    %
    % History (in reverse chronological order)
    % Version:   2.0
    % Author:    Graham Rogers
    % Date:      March 1994
    % Purpose:   eliminated do loops to improve speed
    % Version:   1.0
    % Author:    Kwok W. Cheung, Joe H. Chow
    % Date:      March 1991
    %
    % ***********************************************************
    jay = sqrt(-1);
    exp_ang = exp(jay*ang);
    % Voltage rectangular coordinates
    V_rect = V.*exp_ang;
    CV_rect=conj(V_rect);
    Y_con = conj(Y);
    %vector of conjugate currents
    i_c=Y_con*CV_rect;
    % complex power vector
    S=V_rect.*i_c;
    S=sparse(diag(S));
    Vdia=sparse(diag(V_rect));
    CVdia=conj(Vdia);
    Vmag=sparse(diag(V));
    S1=Vdia*Y_con*CVdia;
    t1=((S+S1)/Vmag)*volt_red';
    t2=(S-S1)*ang_red';
    J11=-ang_red*imag(t2);
    J12=ang_red*real(t1);
    J21=volt_red*real(t2);
    J22=volt_red*imag(t1);
    if nargout > 3
       Jac11 = J11; clear J11
       Jac12 = J12; clear J12
       Jac21 = J21; clear J21
       Jac22 = J22; clear J22
    else
       Jac11 = [J11 J12;
               J21 J22];
    end
 

function     [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line)
    % Syntax:    [Y,nSW,nPV,nPQ,SB] = y_sparse(bus,line) 
    %
    % Purpose:   build sparse admittance matrix Y from the line data
    %
    % Input:     bus  - bus data
    %            line - line data
    %
    % Output:    Y    - admittance matrix
    %            nSW  - total number of swing buses
    %            nPV  - total number generator buses
    %            nPQ  - total number of load buses
    %            SB   - internal bus numbers of swing bus
    %
    % See also:  
    %
    % Calls:    
    %
    % Called By:   loadflow, form_j, calc

    % (c) Copyright 1994-1996 Joe Chow - All Rights Reserved
    %
    % History (in reverse chronological order)
    %
    % Version:   2.0
    % Author:    Graham Rogers
    % Date:      April 1994
    % Version:   1.0
    % Author:    Kwok W. Cheung, Joe H. Chow
    % Date:      March 1991
    %
    % ************************************************************
    global bus_int

    jay = sqrt(-1);
    swing_bus = 1;
    gen_bus = 2;
    load_bus = 3;

    nline = length(line(:,1));     % number of lines
    nbus = length(bus(:,1));     % number of buses
    r=zeros(nline,1);
    rx=zeros(nline,1);
    chrg=zeros(nline,1);
    z=zeros(nline,1);
    y=zeros(nline,1);

    Y = sparse(1,1,0,nbus,nbus);

    % set up internal bus numbers for second indexing of buses
    busmax = max(bus(:,1));
    bus_int = zeros(busmax,1);
    ibus = [1:1:nbus]';
    for i = 1:nbus
      bus_int(round(bus(i,1))) = i;
    end

    % process line data and build admittance matrix Y
      r = line(:,3);
      rx = line(:,4);
      chrg =jay*sparse(diag( 0.5*line(:,5)));
      z = r + jay*rx;     % line impedance
      y = sparse(diag(ones(nline,1)./z));

    % determine connection matrices including tap changers and phase shifters
      from_bus = round(line(:,1));
      from_int = bus_int(from_bus);
      to_bus = round(line(:,2));
      to_int = bus_int(to_bus);
      tap_index = find(abs(line(:,6))>0);
      tap=ones(nline,1);
      tap(tap_index)=1. ./line(tap_index,6);
      phase_shift = line(:,7);
      tap = tap.*exp(-jay*phase_shift*pi/180);

      % sparse matrix formulation
      iline = [1:1:nline]';
      C_from = sparse(from_int,iline,tap,nbus,nline,nline);
      C_to = sparse(to_int,iline,ones(nline,1),nbus,nline,nline);
      C_line = C_from - C_to;

    % form Y matrix from primative line ys and connection matrices
      Y=C_from*chrg*C_from' + C_to*chrg*C_to' ;
      Y = Y + C_line*y*C_line';
      Gb = bus(:,8);     % bus conductance
      Bb = bus(:,9);     % bus susceptance

    % add diagonal shunt admittances
      Y = Y + sparse(ibus,ibus,Gb+jay*Bb,nbus,nbus);

    if nargout > 1
      % count buses of different types
      nSW = 0;
      nPV = 0;
      nPQ = 0;
      bus_type=round(bus(:,10));
      load_index=find(bus_type==3);
      gen_index=find(bus_type==2);
      SB=find(bus_type==1);
      nSW=length(SB);
      nPV=length(gen_index);
      nPQ=length(load_index);
    end

    return


