clear variables;
close all;

% User Options
% fname = 'bdb_v1_4'; % file name to save database
fname = 'test5'; % file name to save database

nb = 3;  % Number of batches to generate
nb_PRBS = 0;  % Number of PRBS batches

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
samt = .5; % sampling time [h]

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
Kcw = 0.05;
% Kcw = 5.5;

% Calculated paramters
nFE= Tf/samt;

% Generate Database Object
bdb = BatchDatabase(size(X0,2), 6, 4, 3, Tf, samt);
bdb.AddICStd(X_std');
bdb.AddMeasStd(Noise_std');

% temporary variables for storing batch data
X = zeros(nFE+1,bdb.nx); % allocate memory for states
U = zeros(nFE, bdb.nu); % allocate memory for inputs
Y = zeros(nFE+1, bdb.ny);
Q = zeros(nFE+1, bdb.nq);

% NOTE: this is considered an output for the process but is controlled by a
% PI controller based on reaching the batch temp setpoint.
jacket_temp_ideal = zeros(nFE,1);
jacket_temp_implemented = zeros(nFE,1);

Terror = zeros(nFE,1);



ODE_inputs = [X0(2) 
    Nom_feed_rate % Set in loop (on when below threshold)
    Aeration_rate
    Agitator_power
    Feed_temp
    Cooling_feed
    300   % this will be set by the PI controller 
    pH];

% %Parameters for PRBS
Tsw = 6;
seed = 7;

%Minimum prbs signal
u_prbs_min = [-0.01
    -0.01
    -0.05
    -1
    -3
    -0.01
    -2
    -0.3]*0.1;
% u_prbs_min(1) Dissolved O2
% u_prbs_min(2) Nom_feed_rate
% u_prbs_min(3) Aeration Rate
% u_prbs_min(4) Agitator_power
% u_prbs_min(5) Feed_temp
% u_prbs_min(6) Cooling_feed
% u_prbs_min(7) Temp jacket. This will be set by the PI controller 
% u_prbs_min(8) pH

%Maximum prbs signal
u_prbs_max = [0.01
    0.01
    0.05
    1
    3
    0.01
    2
    0.3]*0.1;
% u_prbs_max(1) Dissolved O2
% u_prbs_max(2) Nom_feed_rate
% u_prbs_max(3) Aeration Rate
% u_prbs_max(4) Agitator_power
% u_prbs_max(5) Feed_temp
% u_prbs_max(6) Cooling_feed
% u_prbs_max(7) Temp jacket. This will be set by the PI controller 
% u_prbs_max(8) pH


%Create temporary vector to allocate prbs signals
PRBS_signal = zeros(Tf/samt+1,size(ODE_inputs,2));

% Generate prbs signal for each input at a time
for i=1:size(ODE_inputs,1)

    PRBS_signal(:,i) =  gbngen(Tf/samt+1,Tsw,seed,u_prbs_min(i),u_prbs_max(i));
end

% Generate Nominal Batch
feed_started = false;
opts = odeset('reltol', 1e-9);

% Initialize batch
X(1,:) = X0;
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
  
  meas_noise = bdb.GenNoise(1342)';  % Add measurement noise
for k = 1:nFE
    
    % USE PI CONTROL TO DETERMINE THE JACKET TEMPERATURE
    Terror(k) = (Y(k,2)-T_setpt);
    if(k == 1)
        jacket_temp_ideal(k) = Kc*Terror(k) + 298;
    else
        delta_u = Kc*(Terror(k) - Terror(k-1) + samt/tauI*Terror(k));
        jacket_temp_ideal(k) = jacket_temp_ideal(k-1) + delta_u;
    end
    jacket_temp_implemented(k) = jacket_temp_ideal(k);
    
    if(jacket_temp_ideal(k)<290)
       jacket_temp_implemented(k) = 290; 
    elseif(jacket_temp_ideal(k)>323)
        jacket_temp_implemented(k) = 323;
    end
    ODE_inputs(7) =  jacket_temp_implemented(k);
%       ODE_inputs(7) = 298;
    
    % DETERMINE SUBSTRATE FEED RATE ( turn on when below minimum substarte
    % concentration
    if(~feed_started && X(k,1)<Substrate_switch_conc)
        feed_started = true;
    end
    if(feed_started)
       ODE_inputs(2) =  Nom_feed_rate;
    else
        ODE_inputs(2) = 0;
    end
    
    
    % INTEGRATE SYSTEM
    [T_out,X_out]=ode15s('pengsimv2_1',[(k-1) k]*samt,X(k,:),opts,ODE_inputs); 
    
    % STORE RESUTLS IN TEMPORARY MATRICIES
    X(k+1,:) = X_out(end,:);
    
    Y(k+1,:) = [jacket_temp_implemented(k)
        X_out(end,7)
        X_out(end,2)
        X_out(end,6)
        X_out(end,5)
        ODE_inputs(5)] + meas_noise(k,:)';
  
    
    Q(k+1,:) = [X_out(end,1)
         X_out(end,3)
          X_out(end,4)];
      
    U(k,:) = [ODE_inputs(2)
        ODE_inputs(3)
        ODE_inputs(4)
        ODE_inputs(6)];    
end

% % Include the Heat in the outputs
% Y(:,7)=UA*(Y(:,1)-Y(:,2));

xnorm=X;
ynorm=Y;
unorm=U;
qnorm=Q;



% save nominal batch to database
bdb.AddNomBatch(X', Y', U', Q');

figure;
subplot 321
plot(bdb.data.ynom(1,:));
title('Y1');

subplot 322
plot(bdb.data.ynom(2,:));
title('Y2');

subplot 323
plot(bdb.data.ynom(3,:));
title('Y3');

subplot 324
plot(bdb.data.ynom(4,:));
title('Y4');

subplot 325
plot(bdb.data.ynom(5,:));
title('Y5');

subplot 326
plot(bdb.data.ynom(6,:));
title('Y6');

fprintf(['--> Finished with the nominal batch\n']);


% Genearate Initial Conditions for remaining batches
ICs = bdb.GetICs(nb, 7);
ICs = abs(ICs);

figure;
for b = 1:nb  % for each batch
    
    % ZERO VARIABLES BEFORE USING THEM!
    X = zeros(nFE+1,bdb.nx); % allocate memory for states
    U = zeros(nFE, bdb.nu); % allocate memory for inputs
    Y = zeros(nFE+1, bdb.ny);
    Q = zeros(nFE+1, bdb.nq);
    Terror = zeros(nFE,1);
    Verror = zeros(nFE,1);
    Nom_feed_rate1 = zeros(nFE,1);
    O2error = zeros(nFE,1);
    Aeration_rate1 = zeros(nFE,1);
    CO2error = zeros(nFE,1);
    Agitator_power1 = zeros(nFE,1);
    
    
    % initialize states
    X(1,:) = ICs(b,:);
    Y(1,:) = [X(1,7)
        X(1,7)
        X(1,2)
        X(1,6)
        X(1,5)
        ODE_inputs(5)];

    Q(1,:) = [X(1,1)
        X(1,3)
        X(1,4)];
    
    meas_noise = bdb.GenNoise(b*7777)';  % Add measurement noise
    for k = 1:nFE % for each sampling instant

          % USE PI CONTROL TO DETERMINE THE JACKET TEMPERATURE
        Terror(k) = (Y(k,2)-T_setpt);
        if(k == 1)
            jacket_temp_ideal(k) = Kc*Terror(k) + 298;
        else
            delta_u = Kc*(Terror(k) - Terror(k-1) + samt/tauI*Terror(k));
            jacket_temp_ideal(k) = jacket_temp_ideal(k-1) + delta_u;
        end
        jacket_temp_implemented(k) = jacket_temp_ideal(k);

        if(jacket_temp_ideal(k)<290)
           jacket_temp_implemented(k) = 290; 
        elseif(jacket_temp_ideal(k)>323)
            jacket_temp_implemented(k) = 323;
        end
        ODE_inputs(7) =  jacket_temp_implemented(k);
        ODE_inputs(6) = Kcw*(Y(k,2)-ynorm(k,2))+0.06;
        
        ODE_inputs(7) = 309;
%         ODE_inputs(6) = 0.005;
        
        % Use PI control (u1) to determine fementer volume (Y5)
        Kv=-0.01;
        tauI1=1000;
        Verror(k)=Y(k,5)-ynorm(k,5);
        if (k==1)
            Nom_feed_rate1(k)=Kv*Verror(k)+Nom_feed_rate;
        else
            delta_u1 = Kv*(Verror(k) - Verror(k-1) + samt/tauI1*Verror(k));
            Nom_feed_rate1(k)=Nom_feed_rate1(k-1)+delta_u1;
        end
            
            

        % DETERMINE SUBSTRATE FEED RATE ( turn on when below minimum substarte
        % concentration
        if(~feed_started && X(k,1)<Substrate_switch_conc)
            feed_started = true;
        end
        
        if(feed_started)
           ODE_inputs(2) =  Nom_feed_rate1(k);
        else
            ODE_inputs(2) = 0;
        end
        
        
        % Use PI control (u2) to determine dissolved O2 (Y3)
        Kar=-0.4;
        tauI2=1000;
        O2error(k)=Y(k,3)-ynorm(k,3);
        if (k==1)
            Aeration_rate1(k)=Kar*O2error(k)+Aeration_rate;
        else
            delta_u2 = Kar*(O2error(k) - O2error(k-1) + samt/tauI2*O2error(k));
            Aeration_rate1(k)=Aeration_rate1(k-1)+delta_u2;
        end
        ODE_inputs(3)=Aeration_rate1(k);
        
        % Use PI control (u3) to determine dissolved CO2 (Y4)
        Kap=-0.4;
        tauI3=1000;
        CO2error(k)=Y(k,4)-ynorm(k,4);
        if (k==1)
            Agitator_power1(k)=Kap*CO2error(k)+Agitator_power;
        else
            delta_u3 = Kap*(CO2error(k) - CO2error(k-1) + samt/tauI3*CO2error(k));
            Agitator_power1(k)=Agitator_power1(k-1)+delta_u3;
        end
        ODE_inputs(4)=Agitator_power1(k);
        

        % INTEGRATE SYSTEM
        [T_out,X_out]=ode15s('pengsimv2_1',[(k-1) k]*samt,X(k,:),opts,ODE_inputs); 

        % STORE RESUTLS IN TEMPORARY MATRICIES
        X(k+1,:) = X_out(end,:);

        Y(k+1,:) = [jacket_temp_implemented(k)
        X_out(end,7)
        X_out(end,2)
        X_out(end,6)
        X_out(end,5)
        ODE_inputs(5)] + meas_noise(k,:)';

        Q(k+1,:) = [X_out(end,1)
             X_out(end,3)
              X_out(end,4)];

        U(k,:) = [ODE_inputs(2)
            ODE_inputs(3)
            ODE_inputs(4)
            ODE_inputs(6)];    
    end
    
%     % Include the Heat in the outputs
%     Y(:,7)=UA*(Y(:,1)-Y(:,2));


subplot 321
plot(Y(:,1));
title('Y1');
hold on;

subplot 322
plot(Y(:,2));
title('Y2');
hold on;

subplot 323
plot(Y(:,3));
title('Y3');
hold on;

subplot 324
plot(Y(:,4));
title('Y4');
hold on;

subplot 325
plot(Y(:,5));
title('Y5');
hold on;

subplot 326
plot(Y(:,6));
title('Y6');
hold on;
  
    bdb.AddBatch(X', Y', U', Q');
    
    fprintf(['-- Finished batch ' num2str(b) ' of ' num2str(nb) '\n']);
end

%Get batches with PRBS
% Genearate Initial Conditions for PRBS batches
ICs = bdb.GetICs(nb_PRBS, 12);
ICs = abs(ICs);

for b = 1:nb_PRBS  % for each batch
    
    % ZERO VARIABLES BEFORE USING THEM!
    X = zeros(nFE+1,bdb.nx); % allocate memory for states
    U = zeros(nFE, bdb.nu); % allocate memory for inputs
    Y = zeros(nFE+1, bdb.ny);
    Q = zeros(nFE+1, bdb.nq);
    Terror = zeros(nFE,1);

    
    % initialize states
    X(1,:) = ICs(b,:);
    Y(1,:) = [X(1,7)
        X(1,7)
        X(1,2)
        X(1,6)
        X(1,5)
        ODE_inputs(5)];

    Q(1,:) = [X(1,1)
        X(1,3)
        X(1,4)];
    
    meas_noise = bdb.GenNoise(b*7798)';  % Add measurement noise
    
    for k = 1:nFE % for each sampling instant

          % USE PI CONTROL TO DETERMINE THE JACKET TEMPERATURE
        Terror(k) = (Y(k,2)-T_setpt);
        if(k == 1)
            jacket_temp_ideal(k) = Kc*Terror(k) + 298;
        else
            delta_u = Kc*(Terror(k) - Terror(k-1) + samt/tauI*Terror(k));
            jacket_temp_ideal(k) = jacket_temp_ideal(k-1) + delta_u;
        end
        
        % NOTE WE ALSO ADD THE PRBS SIGNAL HERE SO THAT IT IS DONE BEFORE
        % THE CONSTRAINT CHECK.
        jacket_temp_implemented(k) = jacket_temp_ideal(k)+ PRBS_signal(k,7); % Signal for Temp jacket;

        if(jacket_temp_ideal(k)<290)
           jacket_temp_implemented(k) = 290; 
        elseif(jacket_temp_ideal(k)>323)
            jacket_temp_implemented(k) = 323;
        end
        %ODE_inputs(7) =  jacket_temp_implemented(k);

        % DETERMINE SUBSTRATE FEED RATE ( turn on when below minimum substarte
        % concentration
        if(~feed_started && X(k,1)<Substrate_switch_conc)
            feed_started = true;
        end
        if(feed_started)
           ODE_inputs(2) =  Nom_feed_rate;
        else
            ODE_inputs(2) = 0;
        end

        % Add the PRBS on top of the original signal.
        ODE_inputs(1) = X(k,2) + PRBS_signal(k,1); % Signal for Dissolved O2
        ODE_inputs(2) = Nom_feed_rate + PRBS_signal(k,2); % Signal for Nom_feed_rate
        ODE_inputs(3) = Aeration_rate + PRBS_signal(k,3); % Signal for Aeration Rate
        ODE_inputs(4) = Agitator_power + PRBS_signal(k,4); % Signal for Agitator_power
        ODE_inputs(5) = Feed_temp + PRBS_signal(k,5); % Signal for Feed_temp
        ODE_inputs(6) = Cooling_feed + PRBS_signal(k,6); % Signal for Cooling_feed
        ODE_inputs(7) = jacket_temp_implemented(k);
        ODE_inputs(8) = pH + PRBS_signal(k,8); % Signal for pH
        
        %Check if one of the signals is less than zero
        for i=1:size(ODE_inputs,2)
            
            if ODE_inputs(i)<=0
                ODE_inputs(i)=0;
            end
            
        end
        
        % INTEGRATE SYSTEM
        [T_out,X_out]=ode15s('pengsimv2_1',[(k-1) k]*samt,X(k,:),opts,ODE_inputs); 

        % STORE RESUTLS IN TEMPORARY MATRICIES
        X(k+1,:) = X_out(end,:);

        Y(k+1,:) = [jacket_temp_implemented(k) % Jacket temperature
            X_out(end,7)    % Reactor temperature
            X_out(end,2)    % Dissolved O2
            X_out(end,6)    % Dissolved CO2
            X_out(end,5)    % fermentor volume
            ODE_inputs(5)] + meas_noise(k,:)'; % Feed temp
    
        Q(k+1,:) = [X_out(end,1)
             X_out(end,3)
              X_out(end,4)];

        U(k,:) = [ODE_inputs(2)
            ODE_inputs(3)
            ODE_inputs(4)
            ODE_inputs(6)];    
    end
    
%     % Include the Heat in the outputs
%     Y(:,7)=UA*(Y(:,1)-Y(:,2));

    bdb.AddPRBSBatch(X', Y', U', Q', Tsw);
    
    fprintf(['-- Finished batch_prbs ' num2str(b) ' of ' num2str(nb_PRBS) '\n']);
end

save(fname, 'bdb');



