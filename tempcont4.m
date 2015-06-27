% Temperature controller for Penicillin Production Simulator, PenSim (pensim2.m)
% PID digital controller in velocity form
%
% I/O structure : Fc=tempcont4(inpt,err,xtc,flh,flc);
%
% Fc : Flow rate of cooling or heating water
% err : vector of error history
% xtc : vector of Temperature history
% flh : vector of Heating water history
% flc : vector of Cooling water history
% inpt : (1-by-9) vector of constants
% 
% Cenk Undey (c) January 18,2000
% Department of Chemical & Env. Engineering	
% Illinois Institute of Technology, Chicago - IL/USA
%
% All rights reserved. Copyright (c) IIT 2002.

function Fc=tempcont4(inpt,err,xtc,flh,flc)

q=inpt(1);%Sampling counter
error2=inpt(2);%Current value of the error
Kcc=inpt(3);%Controller gain for cooling system
ti_c=inpt(4);%Integral const. for cooling system
td_c=inpt(5);%Derivative const. for cooling system
Kch=inpt(6);%Controller gain for heating system
ti_h=inpt(7);%Integral const. for heating system
td_h=inpt(8);%Derivative const. for heating system
b=inpt(9);%flag for heating/cooling switching (b=0 : Cooling, b=1 : Heating)

%------PID Control-----%
if b==0 %Cooling
   if q==1
      Fcc=0.001;
      xtcn1=xtc(q);
      xtcn2=xtc(q);
      errn1=error2;
      flown1=Fcc;
   elseif q==2
      flown1=flc(q-1);
      xtcn1=xtc(q-1);
      xtcn2=xtc(q-1);
      errn1=err(q-1);
   elseif q>=3
      flown1=flc(q-1);
      xtcn1=xtc(q-1);
      xtcn2=xtc(q-2);
      errn1=err(q-1);
   end
elseif b==1 %Heating
   if q==1
      Fch=0.001;
      xtcn1=xtc(q);
      xtcn2=xtc(q);
      errn1=error2;
      flown2=Fch;
   elseif q==2
      flown2=flh(q-1);
      xtcn1=xtc(q-1);
      xtcn2=xtc(q-1);
      errn1=err(q-1);
   elseif q>=3
      flown2=flh(q-1);
      xtcn1=xtc(q-1);
      xtcn2=xtc(q-2);
      errn1=err(q-1);
   end
end

if b==1 %Heating
   Ph=(error2-errn1); %Proportional
   Ih=(error2/ti_h); %Integral
   Dh=[td_h*(xtc(q)-2*xtcn1+xtcn2)]; %Derivative
   deltah=Kch*[Ph+Ih-Dh];%PID
   Fch=flown2+deltah;
   if Fch<=0
      Fch=1e-4;
   end
   Fc=Fch;
elseif b==0 %Cooling
   Pc=(error2-errn1); %Proportional
   Ic=(error2/ti_c); %Integral
   Dc=[td_c*(xtc(q)-2*xtcn1+xtcn2)]; %Derivative
   deltac=Kcc*[Pc+Ic-Dc]; %PID
   Fcc=flown1-deltac;
   if Fcc<=0
      Fcc=1e-4;
   end
   Fc=Fcc;
end
