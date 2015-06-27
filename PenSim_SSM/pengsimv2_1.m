% Penicillin G production simulator ODE set
% This file is called by pensim2.m wrapper file
%
% I/O structure : dy = pengsimv1_6m(t,y,FLAG,inp1);
%
% Inputs:
% t: time span of integration
% y: initial conditions for that time span
% FLAG=[]; Matlab notation
% inp1: (1-by-14) vector of input parameters
%
% Output:
% dy: matrix of time stamps and ODE solutions
%
% Cenk Undey, Gulnur Birol and Ali Cinar (c) Control Group, February 01, 2000
% Illinois Institute of Technology, Chicago - IL/USA
%
% Acknowledgement: We gratefully acknowledge the ideas and insight of Dr. I. Birol on pH equilibrium equations.
%
% All rights reserved. Copyright (c) IIT 2002.


function dy=pengsimv2_0(t,y,FLAG,inp1)
a=inp1(1); %a is Dissolved Oxygen Concentrarion in g/L. In the article, this variable is cL and it is explict in equations (2), (3), (10) and (12).
F=inp1(2); %F is Feed Flow Rate of Substrate in L/h. In the article, this variable is also named F and it is explicit in equations (5), (11), (14) and (17)
fgg=inp1(3); %fgg is Flow Rate of Oxygen. In the article, this variable is named fg and it is explicit in equation (13).
pww=inp1(4); %pww is Power Input. In the article, this variable is named Pw and it is explicit in equation (13).
Tf=inp1(5); % Tf is Feed Temperature of Substrate in K. In the article, this variable is named Tf and it is explicit in equation (17).
Fc=inp1(6); %Fc is Cooling Water Flow Rate in L/h. In the article, this variable is also named Fc and it is explicit in equation (17). 
Tc = inp1(7); % Temperature of the heating/cooling water
pH=inp1(8); %This is pH that will be used to calculate the concentration of the hydrogen ion.

%Proposed Feed system
% %If the Substrate concentration is less than a desired value
% if y(1)<= 2
%     %Then, Feed the system
%     F=inp1(10);
% else
%     F=inp1(2);
% end

sf=600;% Feed Substrate Solution Concentration in g/L that is fed to the fermenter. In the article, this variable is also named Sf and it is explicit in equations (11) and (17).
hydrogenion=10^(-pH); %That is the concentration of the hydrogen ion. This parameter is used in equation (3). 
a1=(1/0.45);% 1/(Yx/s). It is the inverse of the yield constant Yx/s. This term 1/(Yx/s) is in equation (11). 
a2=(1/0.9);% 1/(Yp/s).  It is the inverse of the yield constant Yp/s. This term 1/(Yp/s) is in equation (11).
a3=(1/0.04);% 1/(Yx/o). It is the inverse of the yield constant Yx/o. This term 1/(Yx/o) is in equation (12).
a4=(1/0.2);% 1/(Yp/o). It is the inverse of the yield constant Yp/o. This term 1/(Yp/o) is in equation (12).
a5=(1/7);% Constant relating CO2 to growth in mmol CO2/ g biomass. In the article, this variable is named alpha1 and it is explicit in equation (18).
mx=0.014;% Maintenance Coefficient on Substrate per h. In the article, this variable is also named mx and it is explicit in equation (11).
mo=0.467;% Maintenance Coefficient on Oxygen per h. In the article, this variable is also named mo and it is explicit in equation (12).
mux=0.092;% Maximum specific growth rate per h. In the article, this variable is named ux and it it explicit in equations (2), (3), and (4).
kx=0.15;% Contois Saturation Constant in g/L. In the article, this variable is also named Kx and it is explicit in equations (2) and (3).
mup=0.005;% Specific Rate of Penicillin Production per h. In the article, this variable is named up and it is explicit in equation (10).
kp=0.0002;% Monod saturation constant. In table 2, it is named as Inhibition Constant, Kp. It is explicit in equation (10).
ki=0.1;% Substrate inhibition const. for product formation in g/L. In the article, it is also named Ki and is explicit in equation (10).
mc=4e-7;% Constant relating CO2 to maintenance energy in mmol CO2 / (g biomass.h). In the article, this variable is named alpha2 and is explicit shown in equation (18).
k4=1e-4;% Constant relating CO2 to penicillin production in mmol CO2/(L.h). In the article, this variable is named alpha3 and is explicit in equation (18).
p=3;% It is a constant that is in the exponent in equation (10).
k=0.04;% 1st order decay rate const. for product. Penicilin Hydrolysis rate constant per h. In the article, it is also named K and it is explicit in equation (9).
Ed=50000;% Activation energy for cell death in cal/mol. In the article, it is also named Ed and it is explicit in equation (3) and (7).
kdd=1e33;% Arrhenius constant for cell death. In the article, it is named Kd and explicit in equation (3) and (7).
Eg=5100;% Activation energy for growth in cal/mol. In the article, it is also named Eg and it is explicit in equation (3) and (7).
kdg=7e3;% Arrhenius constant for growth. In the article, it is named Kg and explicit in equation (3) and (7).
ck=(1e-5);% Proportionality constant, present in table 2, in mol H+/g biomass. In the article, it has the symbol gama and is explicit in equation (5). 
k1=1e-10;% Constants for [H+] in mol/L. Constant K1 named in the article. It is explicit in equations (3) and (4).
k2=7e-5;% Constants for [H+] in mol/L. Constant K2 named in the article. It is explicit in equations (3) and (4).
rcp=(1/1500);% 1/(ro*cp) in 1/(cal/L.C) for the bulk in fermenter. In the article, it is named pCp and is explicit in equation (17).
rq1=60;% heat generated by fermentation in cal/g biomass. Yield of heat generation. In the article, it is named rq1 and is explicit in equation (16).
rq2=1.6783e-4;% heat generated by fermentation in cal/g biomass.h. Constant in heat generation. In the article, it is named rq2 and explicit shown in equation (16).
a10=0.6; % Constant b. In the article, it is also named b and explicit shown in equation (17).

% switched dynamics (different heat coefficients for heating and cooling)
if (Tc > y(7))
   a9=8e5; % This variable is Heat Transfer coeff. of cooling/heating liquid in cal/(h.C). In the article, it is named a and explicit shown in equation (17).
else
   a9=1e3;% This variable is Heat Transfer coeff. of cooling/heating liquid in cal/(h.C). In the article, it is named a and explicit shown in equation (17).
end

clstar=a; %Clstar is the Dissolved Oxygen Concentration at saturation of oxygen ib g/L. In table 2 shows the value for cstar=1.16. It is in equation to calculate Kox and Kop (below).

% Oxygen Limitation Decision
kox=2e-2*clstar*(1+tanh(50*(0.58-y(2))))/2; %Kox and Kop are Oxygen Limitation Constants, present in table 2. They are in equations (2), (3).
kop=5e-4*(1+tanh(50*(0.58-y(2))))/2;% Kox and Kop are Oxygen Limitation Constants, present in table 2. They are in equation (10).

%y1=>S; %y2=>CL; %y3=>X; %y4=P %y5=V %y6=CO2 %y7=T
c1=y(1)/(kx*y(3)+y(1)); %This corresponds to the first term of the equation (2).
c2=y(2)/(kox*y(3)+y(2)); %This corresponds to the second term of the equation (2).
c3=y(1)/(kp+y(1)*(1+y(1)/ki)); %This corresponds to the first term of the equation (10).
c4=(y(2)^p)/(kop*y(3)+(y(2)^p)); %This corresponds to the second term of the equation (10).
c6=kdd*(exp(-Ed/(1.987*y(7)))); %This is the fifth term of the equation (3) and the second of the equation (7).
c7=kdg*(exp(-Eg/(1.987*y(7)))); %This is the fourth term of the equation (4) and first of the equation (7).
c8=(1/(1+k1/hydrogenion+hydrogenion/k2)); %This corresponds to the first term of the equations (3) and (4).
c9=y(3)*(mux*c1*c2*c8*c7-F/y(5)-c6); %This is the first term of the equation (5), uX -FX/V.
c5=rq1*c9*y(5)+rq2*y(3)*y(5);%Heat generated by fermentation in calories. This is the equation (16).
c10=F-1e-4*y(5)*(exp(0.05*(y(7)-273))-1);%$$Wrong$$ (poorly described in paper) % Volume correction. This is the equation (14). The last term of this equation is Floss, which is equation (15).
%+d1*Fb*cab/cb+d1*Fa*cab/cb
%Proposed c10: c10=F+d1*Fb*cab/cb+d1*Fa*cab/cb-2.5e-4*y(5)*[exp(0.05*(y(7)-273))-1];
kla=70*sqrt(fgg)*((pww/y(5))^0.4);% fgg=inflow air rate, pww=agitator power. This is the equation (13).
%old B: B=[(1e-14/hydrogenion-hydrogenion)*y(5)-d1*cab*Fa*step1+d1*cab*Fb*step1]/(y(5)+d1*Fb*step1+d1*Fa*step1);%$$Wrong$$ 
%B=((1e-14/hydrogenion-hydrogenion)*y(5)-d1*cab*Fa*step1-d1*cab*Fb*step1)/(y(5)+d1*Fb*step1+d1*Fa*step1);
%%This is the term B in the equation (5). This entire equation corresponds
%%to equation (6). This term is not used anymore.

%---------ODEs-----
dy(1)=-y(3)*(mux*a1*c1*c2*c8*c7+mup*a2*c3*c4+mx-a1*c6)+F*sf/y(5)-y(1)*c10/y(5); %This is the ODE dS/dt, equation (11).

dy(2)=-y(3)*(mux*a3*c1*c2*c8*c7+mup*a4*c3*c4+mo-a3*c6)+kla*(clstar-y(2))-y(2)*F/y(5); %$$Wrong$$ %This is the ODE dCL/dt, equation (12).
%Proposed dy(2): dy(2)=-y(3)*(mux*a3*c1*c2*c8*c7+mup*a4*c3*c4+mo-a3*c6)+kla*(clstar-y(2))-y(2)*c10/y(5);

dy(3)=c9; %$$Wrong$$ %This is de ODE dX/dt, equation (1).
%Proposed dy(3): dy(3)= y(3)*(mux*c1*c2*c8*c7-c6)-y(3)*c10/y(5);

dy(4)=mup*c3*c4*y(3)-k*y(4)-y(4)*F/y(5); %$$Wrong$$ %This is the ODE dP/dt, equation (9).
%Proposed dy(4): dy(4)= mup*c3*c4*y(3)-k*y(4)-y(4)*c10/y(5);

dy(5)=c10; %This is the ODE dV/dt, equation (14).

%dy(6)=a5*y(3)*(mux*c1*c2*c8*c7-F/y(5)-c6+mc)+k4; %$$Wrong$$ %(THIS IS A TYPO: CHANGE THIS FIRST!!!) %This is the ODE dCO2/dt, equation (18).
dy(6)=a5*y(3)*(mux*c1*c2*c8*c7-F/y(5)-c6) +mc*y(3)+k4;
%Proposed dy(6): dy(6)= a5*dy(3)+ mc*y(3)+ k4;

%dy(7)=ck*c9+((-B+sqrt(B^2+4e-14))/2-hydrogenion)/step1;  %This is the ODE
%dH+/dt, equation (5). This equation was excluded from the model. The pH is
%considered constant.

dy(7)=rcp*c5/y(5)+F*(Tf-y(7))/y(5)-(((rcp/y(5))*a9*(Fc^(a10+1))*(y(7)-Tc))/(Fc+(a9*(Fc^a10))/2000)); %This is the ODE dT/dt, equation (17).
%What the paper claims (we think this might be wrong): dy(7)= rcp*y(9)/y(5)+F*(Tf-y(7))/sf-tfl*(((rcp/y(5))*a9*(Fc^(a10+1))*(y(7)-b*Th-(1-b)*Tc))/(Fc+(a9*(Fc^a10))/2000));

%dy(9)=c5;% $$Wrong$$ %This is the equation dQrxn/dt, equation (16).
%Proposed dy(9): dy(9)= rq1*dy(3)*y(5) + rq2*y(3)*y(5);

dy=[dy(1);dy(2);dy(3);dy(4);dy(5);dy(6);dy(7)]; %This put the ODE's in a column vector.
