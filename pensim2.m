% Penicillin G production simulator 
%
% I/O structure : XX = pensim2(Y0,inputs1,step2,tt,phflag,type);
%
% Inputs :
% Y0(1x9) row vector of initial cond's.
% S,cL,X,P,V,CO2,H,T,Q(rxn), respectively.
% example=> Y0=[15 1.16 0.1 0 100 0.5 10^(-5.0) 298 0]; 
% inputs1(1x4) row vector of setpoints for input vars.
% fg1,pw1,F1,Tf0, respectively.
% example=> inputs1=[8 30 0.045 296 5 298];
% step2 : Sampling interval, h
% tt : Simulation time, h
% type=0 => on-off (pH), type=1=> PID (pH), type=[] or not specified => default=PID
%
% Output(s):
% XX (length(tt)x17) matrix.
% 1st column: time stamp
% Columns 2nd thru 5th: inputs (prbs)
% [Aeration rate(fg2), Agitator Power(pw2), Subs. Feed Rate(Fs2), Subs. Feed Temp.(Tf2)];
% Columns 6th thru 14th outputs
% Conc.'s of Subs., DO, Biomass, Pen., Bulk Vol., CO2, pH, Temp., Generated Heat
% Columns 15th thru 17th: manipulated variables.
% fla: Acid flow rate, flb: Base flow rate, flc: Cooling/Heating water flow rate
%
% phflag=0 => Open-loop, phflag=1 => Closed-loop
% tempflag==0 => Open-loop
%
% References:
% 'A Mechanistic Model for Penicillin Fermentation'
% R.K. Bajpai & M. Reuss  
% J.Chem. Tech. Biotechnol., 1980, 30, 333-344.											
%
% 'A Modular Simulation Package for Penicillin Production'
% G. Birol, C. Undey & A. Cinar
% to appear on Computers & Chemical Engng., 2002
%
% Cenk Undey, Gulnur Birol and Ali Cinar (c) Control Group, February 01, 2000
% Modified by Cenk Undey on August 21, 2000.
% Last revision by Cenk Undey on March 01, 2002.
% Department of Chemical & Env. Engineering												  
% Illinois Institute of Technology, Chicago - IL/USA
% Contact info: undey@iit.edu
%
% All rights reserved. Copyright (c) IIT 2002.

function XX = pensim2(Y0,inputs1,step2,tt,tempcont,phcont,faultinfo,type)

t=[0:step2:tt]';

phflag=phcont(1);
tempflag=tempcont(1);

f1=faultinfo(1);
ftype=faultinfo(2);
ft=[faultinfo(3) faultinfo(4)];
magnitude=faultinfo(5);
varno=faultinfo(6);

if nargin < 8
   type =[];
end

if isempty(type),type=1;end

% PRBS Input Sequence Generation
z = prbs(inputs1,tt,step2);
z = z(:,1:tt/step2);

%--Fault Introduction to PRBS---
if f1==1
    fv = varno; % Faulty input signal
    switch ftype
    case 1 % Step change
        for j=1:tt/step2
            if j >= ft(1)/step2 & j <= ft(2)/step2
                z(fv,j)  = z(fv,j) * (1 + magnitude/100);
            else
                z(fv,j)  = z(fv,j);
            end
        end
    case 2 % Drift change
        c=0;
        for j=1:tt/step2
            if j >= ft(1)/step2 & j <= ft(2)/step2
                c = c + step2; % Time constant, 1/h
                z(fv,j)  = z(fv,j) + (0.01 * c * magnitude);
            else
                z(fv,j)  = z(fv,j);
            end
        end
    end
end
%--------------------------

% if f1==1
%    if ftype==1 % Step change
%       if interval(1)==0
%          p1=[];
%       else
%          p1=zeros(1,interval(1)/step2);
%       end
%       p2a=ones(1,(interval(2)-interval(1))/step2);
%       p2=(magnitude/100)*mean(z(varno,:));
%       if interval(2)==tt
%          p3=[p1 p2*p2a];
%          p4=z(varno,:)+p3;
%       else
%          p2b=z(varno,1:(tt-interval(2))/step2);
%          p3=[p1 p2*p2a zeros(1,length(p2b))];
%          p4=z(varno,:)+p3;
%       end
%       z(varno,:)=p4;
%    else % Drift
%       if interval(2)~=tt
%          interval(2)=tt;
%       end
%       p1=zeros(1,interval(1)/step2);
%       tint=1:step2:[interval(2)-interval(1)]+step2;
%       alpha=faultinfo(5)*mean(z(varno,:));
%       p2a=alpha*tint/step2;%ramp change
%       p3=[p1 p2a];
%       p4=z(varno,:)+p3;
%       
%       a1=find(p4<=0.00001);
%       p4(a1)=0;
%       
%       z(varno,:)=p4;
%    end
% end

fg2=z(1,:);pw2=z(2,:);Fs2=z(3,:);Tf2=z(4,:);

%--PRBS of set points for pH and Temp. 
phsp=z(5,:);Tsp=z(6,:);
%--------------------

X=Y0;
a=Y0(2);
X1=[0];%initialization
step1=0.01;%ODE solver step size
Fd=0;%Subs. Feed Initialization
i=0;%Time counter
subs=0.31;%Batch/fed-batch switching threshold
cb=3;%base conc. [Molar]
cab=0;% acid conc. [Molar]
Fb=1e-4;%acid/base solution feed rate L/h
Fa=1e-5;
Fc=eps;
phcont1=type;

if phcont1==1
   cab=cb;
end
flh=eps;%Hot water flow rate initialization
if Y0(8)<(Tsp(1)-1)
   b=1;
else
   b=0;
end
% b => flag for heating/cooling switching 
%						(b=0 : Cooling, b=1 : Heating)

%--Controller settings for temperature PID

% Cooling %
Kcc=tempcont(2);
ti_c=tempcont(3);
td_c=tempcont(4);

% Heating %
Kch=tempcont(5);
ti_h=tempcont(6);
td_h=tempcont(7);

%--Controller settings for Acid/base pH PID

% Acid Flow %
Kca=phcont(2);
ti_a=phcont(3);
td_a=phcont(4);

% Base %
Kcb=phcont(5);
ti_b=phcont(6);
td_b=phcont(7);
d1=0;

for q=1:(tt/step2)
   % The following routine (if-end) is for batch/fed-batch switching
   if subs>=0.3 % Subs. conc. (g/L) threshold for switching
      Fs2(q)=0;
      F=Fs2(q);
   else
      F=Fs2(q);
   end
   
   Fd=[Fd;F];F=Fd(q);fgg=fg2(q);pww=pw2(q);Tf=Tf2(q);
   
   %-----pH control-----%
   phc(q)=-log(Y0(7))/log(10);%pH value history
   PH=phc(q);%current pH value
   Er(q)=phsp(q)-PH;% error history
   error4=Er(q);%current value of the error
   
   if phcont1==0
      cab=cb;
      if error4>=0.01
         d1=1;Fb=1.2*1.2e-4;% base flow rate
         Fa=0;% acid flow rate
      elseif error4<=-0.01
         d1=1;Fa=1e-5;%d1=0;
         Fb=0;
      elseif error4==0
         d1=0;Fb=0;Fa=0;
      end
      flb(q)=Fb;% base flow rate history
      fla(q)=Fa;% acid flow rate history
   elseif phcont1==1 & q==1%Initialization
      Fb=0;Fa=0;%acid/base flow rate initializations
      flb(q)=Fb;% base flow rate history
      fla(q)=Fa;% acid flow rate history
   end
   
   % If cab==0 then no control on pH
   if phflag==0
      cab=0;
   end
   if phcont1==1
       inpph=[q error4 Kca ti_a td_a Kcb ti_b td_b];
       %Fa=0;%%%%%%%to cancel acid addition 
       flb(q)=Fb;%Base flow rate history
       fla(q)=Fa;%Acid flow rate history
       d1=1;%Constant used in ODE (pengsimv1_6.m) acid/base equilibrium
       % If d1=1 PID routine works but it's designed to be
       % used mainly for on-off controller
   end
   %------------------
   
   %----Temperature Control----%
   xtc(q)=Y0(8);%Temperature history
   xtc2=xtc(q);%current temperature
   err(q)=Tsp(q)-xtc2;%error history 
   error2=err(q);%current value of the error
   
   if error2>=1.0
      b=1;
      flh(q)=Fc;%hot water flow rate history
   else
      b=0;
      flc(q)=Fc;%cold water flow rate history
   end
   
   inpt=[q error2 Kcc ti_c td_c Kch ti_h td_h b];
   
   if tempflag==0
      tfl=0;
   else
      tfl=1;
   end
   %------------------

   %----Inputs to the ODE file (pengsimv1_6.m)------%
   %inp1=[a F fgg pww Tf Fc b Fb Fa cab cb step1 d1 tfl 5.0];
   inp1=[a F fgg pww Tf Fc b Fb Fa cab cb step1 d1 tfl];
   
   %Y0 = [Y0(1:6) Y0(8)];
   
   
   %-ODE file call for MAtlab's stiff ODE solver function ode23t
   [T1,Y1]=ode23t('pengsimv1_6m',[t(q):step1:t(q+1)],Y0,[],inp1);

   [m1 n1]=size(Y1);
   [m2 n2]=size(T1);
   Y0=Y1(m1,:);
   X=[X;Y1(m1,:)];X1=[X1;T1(m2,:)];
   subs=Y0(1);%Current substrate conc.
   
   if tfl~=0 % Open loop for temperature
      %Temperature control with PID controller
      Fc=tempcont4(inpt,err,xtc,flh,flc);% Coolant/hot water flow rate
   end
   
   if phcont1==1
      [Fb,Fa]=phcontv3(inpph,Er,phc,flb,fla);% Fb:Base, Fa:Acid flow rate
   end
end

time=X1;
X(:,7)=-log(X(:,7))/log(10); %Conversion to pH from H+ conc.

noise2=randn(length(X),1)*0.002;% White noise to DO
noise6=randn(length(X),1)*0.06;% White noise to CO2
X(:,2)=X(:,2)+noise2;
X(:,6)=X(:,6)+noise6;

% History of the manipulated vars.
if phflag~=0
   fla=fla(1:tt/step2)';
   flb=flb(1:tt/step2)';
   flc=flc(1:tt/step2)';
else
   fla=[];flb=[];flc=[];
end
Fd=Fd(1:tt/step2)';

XX1=[fg2; pw2; Fd; Tf2]';% Matrix of input variables
X=X(1:tt/step2,:);

XX=[time(1:tt/step2) XX1 X fla flb flc];% Matrix of input and output variables
XX(:,14)=XX(:,14)/1000; % Convert Qrxn (cal/h) to kcal/h.
XX=XX+eps;