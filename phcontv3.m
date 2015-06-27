% pH controller for Penicillin Production Simulator, PenSim (pensim2.m)
% PID digital controller in velocity form
%
% I/O structure : [Fb,Fa]=phcontv3(inpph,Er,phc,flb,fla);
%
% Fb, Fa : Flow rates of base and acid solutions, respectively.
% Er : vector of error history
% phc : vector of pH history
% flb : vector of Base flow rate history
% fla : vector of Acid flow rate history
% inpph : (1-by-8) vector of constants
% 
% Cenk Undey (c) January 18,2000
% Department of Chemical & Env. Engineering	
% Illinois Institute of Technology, Chicago - IL/USA
%
% All rights reserved. Copyright (c) IIT 2002.

function [Fb,Fa]=phcontv3(inpph,Er,phc,flb,fla)

q=inpph(1);%Sampling counter
error4=inpph(2);%Current value of the error
Kca=inpph(3);%Controller gain for acid addition system
ti_a=inpph(4);%Integral const. for acid addition system
td_a=inpph(5);%Derivative const. for acid addition system
Kcb=inpph(6);%Controller gain for base addition system
ti_b=inpph(7);%Integral const. for base addition system
td_b=inpph(8);%Derivative const. for base addition system

%----PID Control for pH - ACID/BASE----%

if q==1
   Fbc=0;Fac=0;
   xtpn1=phc(q);
   xtpn2=phc(q);
   errp1=error4;
   flowp1=Fbc;flowp2=Fac;
elseif q==2
   flowp1=flb(q-1);flowp2=fla(q-1);
   xtpn1=phc(q-1);
   xtpn2=phc(q-1);
   errp1=Er(q-1);
elseif q>=3
   flowp1=flb(q-1);flowp2=fla(q-1);
   xtpn1=phc(q-1);
   xtpn2=phc(q-2);
   errp1=Er(q-1);
end

Pph=(error4-errp1);%Proportional
Ipha=(error4/ti_a);%Integral
Iphb=(error4/ti_b);%Integral
Dpha=[td_a*(phc(q)-2*xtpn1+xtpn2)];%Derivative
Dphb=[td_b*(phc(q)-2*xtpn1+xtpn2)];%Derivative

deltapha=Kca*[Pph+Ipha-Dpha];%PID ACID
deltaphb=Kcb*[Pph+Iphb-Dphb];%PID BASE

Fac=flowp2-deltapha;
Fbc=flowp1+deltaphb;
   
if Fbc<0 
   Fbc=0;
elseif Fbc>0.1 %Physical constraint, max(Fbc)=0.1 L/h
   Fbc=0.1;
end

Fb=Fbc;

if Fac<0 
   Fac=0;
elseif Fac>0.01/1000 %Physical constraint, max(Fac)=0.1 L/h
   Fac=0.01/1000;
end

Fa=Fac;

   