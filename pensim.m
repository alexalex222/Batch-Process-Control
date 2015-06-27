% Penicillin G production simulator wrapper file
%
% Output: X (N-by-17) matrix where N = number of measurements (sampling instances)
%     1st column = time stamp (h)
%     Columns 2 through 5 = inputs (prbs)
%        [Aeration rate(fg2), Agitator Power(pw2), Subs. Feed Rate(Fs2), Subs. Feed Temp.(Tf2)];
%     Columns 6 through 14 = ODE solutions
%        Profiles of Subs. conc., DO conc., Biomass conc., Penicillin conc., Culture Vol., CO2 conc., pH, Temp., Generated Heat
%     Columns 15 through 17: manipulated variables
%        fla: Acid flow rate, flb: Base flow rate, flc: Cooling/Heating water flow rate
%
% Cenk Undey, Gulnur Birol and Ali Cinar (c) Control Group, February 01, 2000
% Modified by Cenk Undey on August 21, 2000.
% Last revision by Cenk Undey on March 01, 2002.
% Department of Chemical & Env. Engineering												  
% Illinois Institute of Technology, Chicago - IL/USA
% Contact info: undey@iit.edu
%
% All rights reserved. Copyright (c) IIT 2002.

function X = pensim

Y0 = [15 1.16 0.1 0 100 0.5 10^(-5.0) 298 0];% Nominal values
% Y0(1) Initial substrate concentration
% Y0(2) Initial Dissolved O2
% Y0(3) Initial Biomass concentration
% Y0(4) Initial penicillin concentration
% Y0(5) Initial fermentor volume
% Y0(6) Initial CO2
% Y0(7) Initial hydrogen ion concentration
% Y0(8) Initial fermentor temperature
% Y0(9) Initial Heat generated

inputs1 = [8 30 0.045 296 5 298];% Nominal values
% inputs1(1) Aeration Rate
% inputs1(2) Agitator Power
% inputs1(3) Feed substrate conc
% inputs1(4) Feed temperature
% inputs1(5) Temperature set point
% inputs1(6) pH set point

tempcont = [1 70 0.5 1.6 5 0.8 0.05];% Nominal values
% tempcont(1) Switch for temp control 1=ON, 0=OFF
% tempcont(2) Coolong controller gain
% tempcont(3) Cooling integral time
% tempcont(4) Cooling derivative time
% tempcont(5) Heating controller gain
% tempcont(6) Heating integral time
% tempcont(7) Heating derivative time

phcont = [1 0.0001 8.4 0.125 0.0008 4.2 0.2625];% Nominal values
% phcont(1) Switch for pH control 1=ON, 0=OFF
% phcont(2) Acid controller gain
% phcont(3) Acid integral time
% phcont(4) Acid derivative time
% phcont(5) Base controller gain
% phcont(6) Base integral time
% phcont(7) Base derivative time

faultinfo = [0 2 0 300 -0.001 3];
% e.g., faultinfo = [1 2 0 200 -0.0009 3];%No fault case

% faultinfo(1) Switch for faults 1=ON, 0=OFF (No faults)
% faultinfo(2) Fault type 1=STEP, 2=RAMP
% faultinfo(3) Starting time (h)
% faultinfo(4) Ending time (h)
% faultinfo(5) Fault magnitude (percent +/-, e.g. -20, 0.01)
% faultinfo(6) Variable number that fault to be introduced
% varno=1:Aeration, 2:Agitator, 3:Subs. Feed

tt = 300;%total production hours to be simulated (h)
step2 = 0.5;%sampling interval (h)

X = pensim2(Y0,inputs1,step2,tt,tempcont,phcont,faultinfo);