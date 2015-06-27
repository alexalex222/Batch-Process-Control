clear;
close all;
clc;



% load state space model and quality model
load SSM_Quality_Scld;

% load batch data
load test9;


% [quality_result, percent_imp, x, y,q, u, u_all, handle, JMPC, time] = QBMPC_EXAMPLE(bdb, bdb, 1, Q_model, mod);
batch_index=50;

% Take the nominal quality as the target if no target is given
qual_target_final = bdb.data.qnom(:,end)';
% qual_target=Q_model.y_mean;

% gams output paramters
gamso.show ='invisible';
gamso.form = 'full';
gamso.compress = true;

%%%% pre-allocate trajectory variables %%%%
x = zeros(bdb.nFE+1, bdb.nx);
y = zeros(bdb.nFE+1, bdb.ny);
q = zeros(bdb.nFE+1, bdb.nq);

% u = bdb.data.unom';
u = bdb.data.u(:,:,batch_index)';


%%%% Initialize trajectory variables %%%%
x(:,:) = bdb.data.x(:,:,batch_index)';
y(:,:) = bdb.data.y(:,:,batch_index)';
q(:,:) = bdb.data.q(:,:,batch_index)';
% note u(1,:) is the first control input and must be calculated

%%%% Handle input constraints and phases %%%%

% Input constraints (columns for phase, rows for varia bles):

% umin = [0;    7.9; 29];
% umax = [0.09; 8.1; 31 ];

umin = [-4; -5; -5; -5];
umax = [4; 5;  5; 5];


t=0.5;
x_now= x(1,:);
u_mpc=[0.045; 8; 30; 0.06];
Terror_past=0;
jacket_temp_ideal_past=298;
[y_next, x_next, Terror_current, jacket_temp_ideal_current]=...
    PenicillinSimulator(t,x_now, u_mpc,Terror_past,jacket_temp_ideal_past);

jacket_temp_ideal=zeros(bdb.nFE+1,1);

x_typ_idx = SSM.x_typ_idx;
y_typ_idx = SSM.y_typ_idx;
u_typ_idx = SSM.x_typ_idx;

% Retrieve the informatio from the subspace model
J_k = SSM.J_k;
C_hat = SSM.C;
PastHorizon = SSM.PastHorizon;
Coeff_temp = (SSM.C*SSM.J_k)';
Coeff_ARX = Convert3DCoeff_ARX(C_hat,J_k,PastHorizon);


% J_k = wgdxvar('J_k', J_k, 2);
% C_hat = wgdxvar('C_hat', C_hat, 2);
if PastHorizon == 1
    Coeff_ARX = wgdxvar('Coeff_ARX', Coeff_ARX, 2);
else
    Coeff_ARX = wgdxvar('Coeff_ARX', Coeff_ARX, 3);
end

%%%% Determine Input Constraints %%%%
      
input_const=[umin umax];
    
input_const = wgdxvar('input_const', input_const, 2);

%%%% Pass the parameter of the quality model %%%%
    
% The bias term is very small so that it can be neglected
Coeff_Quality = SSM.Coeff_Quality(1:end,:);
Coeff_Quality_temp = Convert3DCoeff_Quality(Coeff_Quality, J_k, PastHorizon,...
    size(SSM.y_typ_idx,2),size(SSM.u_typ_idx,2));
% if PastHorizon == 1
%     Coeff_Quality = wgdxvar('Coeff_Quality', Coeff_Quality, 2);
% else
%     Coeff_Quality = wgdxvar('Coeff_Quality', Coeff_Quality, 3);
% end
    
q_meas_mean = SSM.mt_final_mean;
q_meas_std = SSM.mt_final_std;
  
q_pred_mean = SSM.q_final_mean;
q_pred_std = SSM.q_final_std;

q_meas_mean = wgdxvar('q_meas_mean', q_meas_mean, 1);
q_meas_std = wgdxvar('q_meas_std', q_meas_std,1);
    
q_pred_mean = wgdxvar('q_pred_mean', q_pred_mean,1);
q_pred_std = wgdxvar('q_pred_std', q_pred_std,1);

Ymean = SSM.Ymean;
Ysigma = SSM.Ysigma;
Umean = SSM.Umean;
Usigma = SSM.Usigma;
    
    
%%%% Get Desired Quality %%%%

qual_target=qual_target_final.*[0.8 1.5 1.2];

if SSM.ScaldQualityModel == 1
    qual_set = (qual_target-SSM.q_final_mean)./SSM.q_final_std;
else
    qual_set = qual_target;
end
qual_set = wgdxvar('qual_set', qual_set, 1);


    

for i = 3:bdb.nFE
    
   
    
    % Determine the prediction length based on current time
    PredictionStep = bdb.nFE-PastHorizon+1;
    
    
    
%     num_initial = size(Q_model.x0_index,2);
%     num_meas = size(Q_model.x_index,2)+size(Q_model.y_index,2) + size(Q_model.u_index,2);
    %u(i:end,:) = bdb.unom(:,i:end)'%%% TEMPORARY REMOVE THIS CODE %%%
    
    %%%% create gams sets %%%%
    
    % time steps to batch termination T
    T.name = 'T';
    for j = 1:PredictionStep+2-i+PastHorizon
        T.uels{j} = j;
    end
    
    % Measured Variables M
    M.name = 'M';
    for j = 1:(size(SSM.y_typ_idx,2)+size(SSM.u_typ_idx,2))
       M.uels{j} = j; 
    end
    
    % Predicted Variables P
    P.name = 'P';
    for j = 1:size(SSM.y_typ_idx,2)
       P.uels{j} = j; 
    end
    
    % lag in the past horizon
    L.name = 'L';
    for j = 1:PastHorizon
        L.uels{j} = j; 
    end

    % Qualities Q
    Q.name = 'Q';
    for j = 1:bdb.nq
        Q.uels{j} = j;
    end
    
    % Dimension of final state
    F.name = 'F';
    for j = 1:size(SSM.A,1)
        F.uels{j} = j;
    end
    
   
    
    % set used only for passing input constraints
    max_min.name = 'max_min';
    max_min.uels{1} = 1;
    max_min.uels{2} = 2;
    
    % Input U (used only for passing input constraints)
    U.name = 'U';
    for j = 1:size(SSM.u_typ_idx,2)
        U.uels{j} = j;
    end
    
    
    

    wgdx('QBMPCv2sets', T, M, P, L,  Q, F, max_min, U);
    
    
    %%%% Pass current data to GAMs %%%%
    % copy in all inputs (from nominal in first time step then from MPC
    data(:,size(SSM.y_typ_idx,2)+1:size(SSM.y_typ_idx,2)+ size(SSM.u_typ_idx,2)) = ...
        u(i-PastHorizon:end,SSM.u_typ_idx);
    
    % copy intial conditons
    data(1:PastHorizon,1:size(SSM.y_typ_idx,2)) = y(i-PastHorizon:i-1,SSM.y_typ_idx);
    
    [d1,~] = size(data);
    Coeff_Quality = zeros(size(SSM.y_typ_idx,2)+ size(SSM.u_typ_idx,2),bdb.nq,d1);
%     Coeff_Quality = 0.00000001*ones(size(SSM.y_typ_idx,2)+ size(SSM.u_typ_idx,2),bdb.nq,d1);
    Coeff_Quality(:,:,end-1:end) = Coeff_Quality_temp;
    
    Coeff_Quality = wgdxvar('Coeff_Quality', Coeff_Quality, 3);
    
    data = ScaleData(data, [Ymean Umean], [Ysigma Usigma]);
    
    data = wgdxvar('data', data, 2);

    
    
    
    
    
   
   
    
    
    wgdx('MPCinitialize', data, Coeff_ARX, ...
            q_meas_mean, q_meas_std, q_pred_mean, q_pred_std, qual_set, input_const, Coeff_Quality);
    
    
    
    
    tic
    gams('QBMPC_final_state nlp=conopt lo=0');
    time(i) = toc;
    
    s1.name = 'JMPC';
    s1.form = 'full';
    s1.compress = true;
    JMPC_out = rgdx('matsol2', s1);
    
    s2.name = 'X';
    s2.form = 'full';
    s2.compress = true;
    x_out = rgdx('matsol2', s2);
    
    s3.name = 'QUALITY';
    s3.form = 'full';
    s3.compress = true;
    q_pred = rgdx('matsol2', s3);
    
    
        fprintf(['*** Batch ' num2str(batch_index) ' ***\n'])
    
    fprintf(['  i = ' num2str(i) ' JMPC = ' num2str(JMPC_out.val) '\n']);
    
    JMPC(i) = JMPC_out.val;
    
    x_out = x_out.val;
    x_out = ScaleBack(x_out, [Ymean Umean], [Ysigma Usigma]);
    q_pred = q_pred.val;
    if size(x_out,2)<size(SSM.y_typ_idx,2)+size(SSM.u_typ_idx,2)
        x_out(:,size(x_out,2)+1:size(SSM.y_typ_idx,2)+size(SSM.u_typ_idx,2))=0;
    end

%     q_pred = q_pred'.*Q_model.y_std + Q_model.y_mean;
    
%     q_pred_traj(i,:) = q_pred(:);
    
    % extract input moves (for the rest of the batch) from the gams output
    u(i:end,:) = x_out(3:end,size(SSM.y_typ_idx,2)+1:end);
    
    
    
   
    
    %%%% Simulate the next step of batch with the given input trajectory
    
    u_mpc=u(i,:);
    [y_next, x_next,q_next, Terror_current, jacket_temp_ideal_current,noise]=...
    PenicillinSimulator(t,x_now, u_mpc,Terror_past,jacket_temp_ideal_past);

    
    x(i+1,:) = x_next;
    y(i+1,:) = y_next;
    q(i+1,:) = q_next;
    
    x_now=x_next;
    jacket_temp_ideal(i,:)=jacket_temp_ideal_current;
    Terror_past=Terror_current;
    jacket_temp_ideal_past=jacket_temp_ideal_current;
    
    %%%% Epty GAMs variables %%%%
    clear data Coeff_Quality;
    
%     clear Coeff_ARX;
    
    %%%% Empty GAMs Sets %%%%
    clear T M P L Q U max_min;
    

end
quality_result = q(end,1:3);
%%
q_old = bdb.data.q(1:3,end,batch_index)';

percent_imp = (quality_result - q_old)./(qual_target - q_old) * 100;
for i = 1:size(percent_imp,2)
   if(percent_imp(i)>100)
       percent_imp(i) = 200 - percent_imp(i);
   end
end