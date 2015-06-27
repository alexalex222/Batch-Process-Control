
function [y_next, x_next, q_next, Terror_current, jacket_temp_ideal_current, meas_noise]=...
    PenicillinSimulator(t,x_now, u_mpc,Terror_past,jacket_temp_ideal_past,flag)



% t=0.5;
% x_now= [15 1.16 0.1 0 100 0.5 298];
% u_mpc=8;
% Terror_past=0;
% jacket_temp_ideal_past=298;

% This function is to calculate the next setp output of the penicllin
% fermentation system given the u form MPC

% Input:
% t: time
% x_now: current system state. The state is provied from the previous step
% of simulation. In reality, the state cannot be measuured.
% u_mpc: the control input from MPC. u1: Glucose feed rate; u2: Aeration rate;
% u3: agitatior power;
% flag: notation.
% Terror_past: past jakect temperature error (from past step).
% jacket_temp_ideal_past: past ideal jacket temperature (from past step)

% Output:
% y_next: measurement in the next step.
% x_next: system state in the next step. Cannot be measured in reality.
% Terror_current: current temperature error (for next step)
% jacket_temp_ideal_current: current ideal jacket temperatue (for next step)


% Standard Deviations in Initial Conditions
X_std = [.01 .07 .005 .001 4 .05 .005];
% (1) Substrate (glucose) concentration
% (2) Dissolved O2
% (3) Biomass concentration
% (4) penicillin concentration
% (5) fermentor volume
% (6) CO2
% (7) fermentor temperature

% Standard Devaitons in Noise
Noise_std = [.01 .01 .001 .003 .01 .01]*0.01;
% (1) Temp Jacket
% (2) Temp Reactor
% (3) O2
% (4) CO2
% (5) Volume
% (6) Glucose Feed Temp

Tf = 300;  % Total batch duration [h]
samt = 0.5; % sampling time [h]




% Simulation steps
nFE= Tf/samt;

% Operating Parameters
Substrate_switch_conc = 0.3; % Subs. conc. [g/L] threshold for switching
Nom_feed_rate= 0.0450; % Feed rate for fed batch part
Agitator_power =30; %Power Input [?]. In the article, this variable is named Pw and it is explicit in equation (13). This value comes from pensim line 34.
Aeration_rate= 8; % Flow Rate of Oxygen. In the article, this variable is named fg and it is explicit in equation (13). This value comes from pensim line 34.
Cooling_feed = .06; %Cooling Water Flow Rate in L/h. In the article, this variable is also named Fc and it is explicit in equation (17). This value comes from pensim line 149.


% Nominal Sates
 X0 = [15 1.16 0.1 0 100 0.5 298]; %Initial conditions for the states.
% X0(1) Initial substrate (glucose) concentration
% X0(2) Initial Dissolved O2
% X0(3) Initial Biomass concentration
% X0(4) Initial penicillin concentration
% X0(5) Initial fermentor volume
% X0(6) Initial CO2
% X0(7) Initial fermentor temperature

% % Heat Transfer Parameter
% UA = 1000; % Heat transffer coefficient [cal/h.C]

% Fixed Inputs:
pH = 5.0;
Feed_temp = 296; % Feed Temperature of Substrate [K]. In the article, this variable is named Tf and it is explicit in equation (17). This value comes from pensim line 34.
T_setpt = 298; % Set point for the batch temperature [K]

% Temperature PI tuning parameters
Kc = -4000;
tauI = 1000;

% Generate Database Object
bdb = BatchDatabase(size(X0,2), 6, 6, 3, t, samt);
bdb.AddICStd(X_std');
bdb.AddMeasStd(Noise_std');





% ODE_inputs = [X0(2) 
%     Nom_feed_rate % Set in loop (on when below threshold)
%     Aeration_rate
%     Agitator_power
%     Feed_temp
%     Cooling_feed
%     300   % this will be set by the PI controller 
%     pH];

ODE_inputs = [X0(2) 
    u_mpc(1) % feed rate Set in loop (on when below threshold)
    u_mpc(2)  %Aeration_rate
    u_mpc(3)
    Feed_temp
    u_mpc(4)
    300   % this will be set by the PI controller 
    pH];


% Generate Nominal Batch
feed_started = false;
opts = odeset('reltol', 1e-9);

% Initialize current batch state
X(1,:) = x_now;
% Y(1,:) = [X(1,7)
%         X(1,7)
%         X(1,2)
%         X(1,6)
%         X(1,5)
%         ODE_inputs(5)];
Y(1,:) = [X(1,7)
        X(1,7)
        X(1,2)
        X(1,6)
        X(1,5)
        ODE_inputs(5)];
    
    
 Q(1,:) = [X(1,1)
     X(1,3)
      X(1,4)];
  
 
   
  meas_noise = randn(6, 1);  % Add measurement noise
  meas_noise = meas_noise.*Noise_std'; 
%     b=randi(100);
%     meas_noise = bdb.GenNoise(b*7777)';
  
    
    % USE PI CONTROL TO DETERMINE THE JACKET TEMPERATURE
    Terror_current = Y(:,2)-T_setpt;
    if(x_now == X0)
        jacket_temp_ideal_current = Kc*Terror_current + 298;
    else
        delta_u = Kc*(Terror_current - Terror_past + samt/tauI*Terror_current);
        jacket_temp_ideal_current = jacket_temp_ideal_past + delta_u;
    end
    jacket_temp_implemented = jacket_temp_ideal_current;
    
    if(jacket_temp_ideal_current<290)
       jacket_temp_implemented = 290; 
    elseif(jacket_temp_ideal_current>323)
        jacket_temp_implemented = 323;
    end
    ODE_inputs(7) =  jacket_temp_implemented;
    
    % DETERMINE SUBSTRATE FEED RATE ( turn on when below minimum substarte
    % concentration
    if(~feed_started && x_now(:,1)<Substrate_switch_conc)
        feed_started = true;
    end
    
    if(feed_started)
       ODE_inputs(2) =  u_mpc(1);
    else
        ODE_inputs(2) = 0;
    end
    
    
    % INTEGRATE SYSTEM
    [~,X_out]=ode15s('pengsimv2_1',t,x_now,opts,ODE_inputs); 
    
    % STORE RESUTLS IN TEMPORARY MATRICIES
    x_next = X_out(end,:);
    
    y_next = [jacket_temp_implemented
        X_out(end,7)
        X_out(end,2)
        X_out(end,6)
        X_out(end,5)
        ODE_inputs(5)] + meas_noise;
  
    
    q_next = [X_out(end,1)
         X_out(end,3)
          X_out(end,4)];
      
    u_next = [ODE_inputs(2)
        ODE_inputs(3)
        ODE_inputs(4)
        ODE_inputs(6)];    
end
  
  


