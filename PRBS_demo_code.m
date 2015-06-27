clear variables;
close all;

Tsw = 6;

seed = 12; % my favorite number

Tf = 300;  % Total batch duration [h]
samt = .5; % sampling time [h]

u_without_PRBS = 1/1500*(1:1:Tf/samt+1).^(1.5); % THIS IS DUMMY CODE TO CREATE A SAMPLE INPUT TRAJECTORY.  You will use whatever the ODE input was originally here.

% You will need to decide what is appropriate for each of the ode inputs
u_prbs_min = [-0.5];
u_prbs_max = [0.5];

% This code generates the PRBS signal
%  (not this creates the signal for only one input at a time.  You will
%  need to do this for each ODE input.
PRBS_signal =  gbngen(Tf/samt+1,Tsw,seed,u_prbs_min,u_prbs_max);

u_to_implement = u_without_PRBS + PRBS_signal;  % Add the PRBS on top of the original signal.

stairs(u_without_PRBS,'r');
hold on;
stairs(u_to_implement);
