clear;
close all;
clc;

load Q_model_Inter1;
load Q_model_Inter2
load Q_model_final;
% load mod;
load MultiMod_inter1;
load MultiMod_inter2;
load MultiMod_final;

load test9;


% [quality_result, percent_imp, x, y,q, u, u_all, handle, JMPC, time] = QBMPC_EXAMPLE(bdb, bdb, 1, Q_model, mod);
batch_index=1;
% Mult_model=mod;
% Take the nominal quality as the target if no target is given
IndexQinter1=201;
IndexQinter2=401;
qual_target_inter1 = bdb.data.qnom(:,IndexQinter1)';
qual_target_inter2 = bdb.data.qnom(:,IndexQinter2)';
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
u = zeros(size( bdb.data.unom'));

u_all = zeros(bdb.nFE, bdb.nu, bdb.nFE);


%%%% Initialize trajectory variables %%%%
x(1,:) = bdb.data.x(:,1,batch_index)';
y(1,:) = bdb.data.y(:,1,batch_index)';
q(1,:) = bdb.data.q(:,1,batch_index)';
% note u(1,:) is the first control input and must be calculated

%%%% Handle input constraints and phases %%%%

% Input constraints (columns for phase, rows for varia bles):

umin = [0;    7.9; 29];
umax = [0.09; 8.1; 31 ];




t=0.5;
x_now= x(1,:);
u_mpc=[0.045; 8; 30];
Terror_past=0;
jacket_temp_ideal_past=298;
[y_next, x_next, Terror_current, jacket_temp_ideal_current]=...
    PenicillinSimulator(t,x_now, u_mpc,Terror_past,jacket_temp_ideal_past);

jacket_temp_ideal=zeros(bdb.nFE+1,1);

for i = 1:bdb.nFE
    
    % determine which local multimodel to use based on current time
    % time
    if i<IndexQinter1
        Mult_model=MultiMod_inter1;
    else if i<IndexQinter2
            Mult_model=MultiMod_inter2;
        else 
            Mult_model=MultiMod_final;
        end
    end
    
    %%%% Retrieve information about the multi-model %%%%
x_typ_idx = Mult_model.getTypeIndex('state');
y_typ_idx = Mult_model.getTypeIndex('output');
u_typ_idx = Mult_model.getTypeIndex('input');


%Mult-model coefficients
if(strcmpi(Mult_model.Regression_Type, 'PLS'))
    PLS_intercept = Mult_model.Coeffs(1,:);
    for c = 1:Mult_model.Opt_Clusters
       
        
       mult_coeffs_intercept(:,c) =  Mult_model.Coeffs((c-1)*(1+Mult_model.NumMeasured)+2,:);
       mult_coeffs(:,:,c) = Mult_model.Coeffs(...
           (c-1)*(1+Mult_model.NumMeasured)+3:c*(1+Mult_model.NumMeasured)+1,:);
    end
else
    PLS_intercept = zeros(1,size(Mult_model.Coeffs,2));
    for c = 1:Mult_model.Opt_Clusters
       mult_coeffs_intercept(:,c) =  Mult_model.Coeffs((c-1)*(1+Mult_model.NumMeasured)+1,:);
       mult_coeffs(:,:,c) = Mult_model.Coeffs(...
           (c-1)*(1+Mult_model.NumMeasured)+2:c*(1+Mult_model.NumMeasured),:);
    end
end

mult_coeffs = wgdxvar('mult_coeffs', mult_coeffs, 3);
mult_coeffs_intercept = wgdxvar('mult_coeffs_intercept', mult_coeffs_intercept, 2);
PLS_intercept = wgdxvar('PLS_intercept', PLS_intercept, 1);

% scaling for mult-model
meas_mean = Mult_model.Measure_Means;
meas_std = Mult_model.Measure_Std;

pred_mean = Mult_model.Predicted_Means;
pred_std = Mult_model.Predicted_Std;


meas_mean = wgdxvar('meas_mean', meas_mean, 1);
meas_std = wgdxvar('meas_std', meas_std, 1);

pred_mean = wgdxvar('pred_mean', pred_mean, 1);
pred_std = wgdxvar('pred_std', pred_std, 1);

% Cluster centers for Mult-Model
Clust_Cent = Mult_model.Centers;

Clust_Cent = wgdxvar('FCMc', Clust_Cent, 2);
    
    % determine which intermedieate quality to predict based on current
    % time
    if i<IndexQinter1
        Q_model=Q_model_Inter1;
    else if i<IndexQinter2
            Q_model=Q_model_Inter2;
        else 
            Q_model=Q_model_final;
        end
    end
    
    % Determine the prediction length based on current time
    if i<IndexQinter1
        PredictionStep=IndexQinter1-1;
    else if i<IndexQinter2
            PredictionStep=IndexQinter2-1;
        else 
            PredictionStep=bdb.nFE;
        end
    end
    
    % Deterimine the quality target based on current time
    if i<IndexQinter1
        qual_target=qual_target_inter1;
    else if i<IndexQinter2
            qual_target=qual_target_inter2;
        else 
            qual_target=qual_target_final.*[0.8 1.3 1.3];
        end
    end
    
    num_initial = size(Q_model.x0_index,2);
    num_meas = size(Q_model.x_index,2)+size(Q_model.y_index,2) + size(Q_model.u_index,2);
    %u(i:end,:) = bdb.unom(:,i:end)'%%% TEMPORARY REMOVE THIS CODE %%%
    
    %%%% create gams sets %%%%
    % time steps to batch termination T
    T.name = 'T';
    for j = 1:PredictionStep+2-i
        T.uels{j} = j;
    end
    
    % Measured Variables M
    M.name = 'M';
    for j = 1:Mult_model.NumMeasured
       M.uels{j} = j; 
    end
    
    % Predicted Variables P
    P.name = 'P';
    for j = 1:Mult_model.NumPredicted
       P.uels{j} = j; 
    end
    
    %Clusters C
    C.name = 'C';
    for j = 1:Mult_model.Opt_Clusters
        C.uels{j} = j;
    end

    % Qualities Q
    Q.name = 'Q';
    for j = 1:size(Q_model.coeffs,2)
        Q.uels{j} = j;
    end
    
    % Number of PCs 
    PC.name = 'PC';
    for j = 1:Q_model.nPC
       PC.uels{j} = j; 
    end
    
    % Past data lenght (used for SPE calc)
    PD.name = 'PD';
    for j = 1:num_initial + num_meas*(i-1)
        PD.uels{j} = j;
    end
    
    % set used only for passing input constraints
    max_min.name = 'max_min';
    max_min.uels{1} = 1;
    max_min.uels{2} = 2;
    
    % Input U (used only for passing input constraints)
    U.name = 'U';
    for j = 1:size(Mult_model.Measurements{u_typ_idx},2)
        U.uels{j} = j;
    end
    
    
    

    wgdx('QBMPCv2sets', T, M, P, C, Q, max_min, U, PC, PD );
    %%%% Pass current data to GAMs %%%%
    % copy in all inputs (from nominal in first time step then from MPC
    data(:,size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)...
        +1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2) ...
        + size(Mult_model.Measurements{u_typ_idx},2)) = ...
        u(i:end,Mult_model.Measurements{u_typ_idx});
    
    % copy intial conditons
    data(1,1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)) = ...
        [x(i,Mult_model.Measurements{x_typ_idx}) y(i,Mult_model.Measurements{y_typ_idx})];
    
    
    % scale everything
    data = data - repmat(Mult_model.Measure_Means, bdb.nFE-(i-1),1);
    data = data./ repmat(Mult_model.Measure_Std, bdb.nFE- (i-1),1);
    
    data1=data;
    data = wgdxvar('data', data, 2);

    %%%% Determine Input Constraints %%%%
    % get intput means and stds
    input_meas_means = Mult_model.Measure_Means(size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)...
        +1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2) ...
        + size(Mult_model.Measurements{u_typ_idx},2));
    input_meas_std = Mult_model.Measure_Std(size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)...
        +1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2) ...
        + size(Mult_model.Measurements{u_typ_idx},2));
    
%     for j = i:bdb.nFE+1
%         phase = find(tphase<=j*bdb.samt, 1,'last');
%         input_const(:,:,j-i+1) = [umin(:,phase) umax(:,phase)];
%         
%         % scale constraint
%         input_const(:,:,j-i+1) = (input_const(:,:,j-i+1) - repmat(input_meas_means',1,2))./repmat(input_meas_std',1,2);
%     end
      
      input_const=[umin umax];
      % Scale constraint
      input_const = (input_const - repmat(input_meas_means',1,2))./repmat(input_meas_std',1,2);
      input_const1= input_const;

    % The dimension is problematic
    input_const = wgdxvar('input_const', input_const, 2);
    
    
    %%%% Determine previous contribution to Quality %%%%
    
    % unfold past data
    past_data = x(1,Q_model.x0_index);
    
    next_entry = size(Q_model.x0_index,2) +1;
    for j = 1:i-1 % NOTE: this goes to i-1 because the current measurements will be handled in GAMs
       past_data(next_entry:next_entry+size(Q_model.x_index,2)-1) = x(j,Q_model.x_index);
       next_entry = next_entry + size(Q_model.x_index,2);
       
       past_data(next_entry:next_entry+size(Q_model.y_index,2)-1) = y(j,Q_model.y_index);
       next_entry = next_entry + size(Q_model.y_index,2);
       
       past_data(next_entry:next_entry+size(Q_model.u_index,2)-1) = u(j,Q_model.u_index);
       next_entry = next_entry + size(Q_model.u_index,2);
    end
    
    % scale past data
    past_data = (past_data - Q_model.x_mean(1:next_entry-1))./Q_model.x_std(1:next_entry -1);
    
    % calculate contribution from past data
    contrib_past = past_data*Q_model.coeffs(1:next_entry -1,:); %+ Q_model.coeffs(1,:);
    
    contrib_past = wgdxvar('Contrib_past', contrib_past, 1);
    
    %%%% Get Quality Model Coeffcients %%%%

    for j = i:PredictionStep+1
       qual_coeffs(:,:,j-i+1) = Q_model.coeffs(num_initial+(j-1)*num_meas+1:num_initial+j*num_meas,:);
    end
    
    qaul_coeffs = wgdxvar('qual_coeffs', qual_coeffs, 3);
    
    %%%% Get Quality Model Scaling %%%%
    for j = i:PredictionStep+1
        q_meas_mean(j-i+1,:) = Q_model.x_mean(num_initial+(j-1)*num_meas+1:num_initial+j*num_meas)';
        q_meas_std(j-i+1,:) = Q_model.x_std(num_initial+(j-1)*num_meas+1:num_initial+j*num_meas)';
    end
    
    q_pred_mean = Q_model.y_mean;
    q_pred_std = Q_model.y_std;

    q_meas_mean = wgdxvar('q_meas_mean', q_meas_mean, 2);
    q_meas_std = wgdxvar('q_meas_std', q_meas_std,2);
    
    q_pred_mean = wgdxvar('q_pred_mean', q_pred_mean,1);
    q_pred_std = wgdxvar('q_pred_std', q_pred_std,1);
    
    %%%% Quality model Tsqrd stuff
    t_past_comp = past_data*Q_model.wstar(1:num_initial + (i-1)*num_meas,:);  % contribution to scores from past measurements
    t_past_comp_gdx = wgdxvar('t_past_comp', t_past_comp,1);
    
    for j = i:PredictionStep+1
        w_star_reformed(:,:,j-i+1) = Q_model.wstar(num_initial+(j-1)*num_meas +1:num_initial + j*num_meas,:);
    end
    w_star_reformed = wgdxvar('w_star', w_star_reformed, 3);
    
     std_scores = Q_model.sa;
    std_scores = wgdxvar('std_scores', std_scores, 1);
    
    Tsqrd_limit = Q_model.T_lim_calc(.95);
    Tsqrd_limit = wgdxvar('Tsqrd_limit', Tsqrd_limit, 0);
    
    %%%% Quality SPE calculation
    past_data_gdx = wgdxvar('past_data', past_data, 1);
    
    w_star_past = Q_model.wstar(1:num_initial+(i-1)*num_meas,:);
    w_star_past = wgdxvar('w_star_past', w_star_past, 2);
    
    SPE_limit = Q_model.SPE_lim_calc(.9999);
    SPE_limit = wgdxvar('SPE_limit', SPE_limit, 0);
    
    %%%% Get Desired Quality %%%%
    qual_set = (qual_target-Q_model.y_mean)./Q_model.y_std;
    qual_set = wgdxvar('qual_set', qual_set, 1);
    
    qual_lo = (0.6 - Q_model.y_mean(1))/Q_model.y_std(1);
    qual_lo = wgdxvar('qual_lo', qual_lo, 0);
   
    %%% get input weights
    % We do not need input weight
%     input_move_weights_wgdx = wgdxvar('input_move_weights', input_move_weights(1,i:end), 1);
    
    
%     wgdx('MPCinitialize', data, mult_coeffs, mult_coeffs_intercept,PLS_intercept, ...
%         Clust_Cent, meas_mean, meas_std, pred_mean, pred_std, contrib_past,...
%         qaul_coeffs, q_meas_mean, q_meas_std, q_pred_mean, q_pred_std, qual_set, input_const, std_scores,...
%         t_past_comp_gdx, w_star_reformed, Tsqrd_limit, past_data_gdx, w_star_past, SPE_limit, qual_lo, input_move_weights_wgdx);
      wgdx('MPCinitialize', data, mult_coeffs, mult_coeffs_intercept,PLS_intercept, ...
        Clust_Cent, meas_mean, meas_std, pred_mean, pred_std, contrib_past,...
        qaul_coeffs, q_meas_mean, q_meas_std, q_pred_mean, q_pred_std, qual_set, input_const, std_scores,...
        t_past_comp_gdx, w_star_reformed, Tsqrd_limit, past_data_gdx, w_star_past, SPE_limit, qual_lo);
    
    
    
    tic
    gams('QBMPC_Pen nlp=conopt lo=0');
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
    q_pred = q_pred.val;
    if size(x_out,2)<7
        x_out(:,size(x_out,2)+1:7)=0;
    end
    x_out = x_out.*repmat(Mult_model.Measure_Std,PredictionStep-i+2,1) + repmat(Mult_model.Measure_Means,PredictionStep-i+2,1);
    q_pred = q_pred'.*Q_model.y_std + Q_model.y_mean;
    
    q_pred_traj(i,:) = q_pred(:);
    
    % extract input moves (for the rest of the batch) from the gams output
    u(i:PredictionStep,1:3) = x_out(1:end-1,5:7);
    
    u_all(i:PredictionStep,:,i) = u(i:PredictionStep,:);
    
    
   
    
    %%%% Simulate the next step of batch with the given input trajectory
    
    u_mpc=u(i,1:3);
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
    clear data contrib_past qaul_coeffs q_meas_mean q_meas_std q_pred_mean q_pred_std qual_set input_const w_star_reformed std_scores Tsqrd_limit t_past_comp;
    
    clear mult_coeffs mult_coeffs_intercept PLS_intercept meas_mean meas_std pred_mean pred_std Clust_Cent;
    
    %%%% Empty GAMs Sets %%%%
    clear T M P C Q U max_min;
    

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