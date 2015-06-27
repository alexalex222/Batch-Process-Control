function [quality_result, percent_imp, x, y,q, u, u_all, handle, JMPC, time] = QBMPC_EXAMPLE(bdb, batch_data, batch_index, Q_model, Mult_model, qual_target)

% Take the nominal quality as the target if no target is given
if(nargin<=5)
    qual_target = bdb.data.qnom(:,end)';
end

% gams output paramters
gamso.show ='invisible';
gamso.form = 'full';
gamso.compress = true;


%%%% pre-allocate trajectory variables %%%%
x = zeros(bdb.nFE+1, bdb.nx);
y = zeros(bdb.nFE+1, bdb.ny);
q = zeros(bdb.nFE+1, bdb.nq);
% u = batch_data.unom';% u = batch_data.unom'; %$$Would it be bdb.data.unom ?
u = bdb.data.unom';

u_all = zeros(bdb.nFE, bdb.nu, bdb.nFE);

%%%% Initialize trajectory variables %%%%
% x(1,:) = batch_data.x(:,1,batch_index);
% y(1,:) = batch_data.y(:,1,batch_index);
% q(1,:) = batch_data.q(:,1,batch_index);
x(1,:) = bdb.data.x(:,1,batch_index)';
y(1,:) = bdb.data.y(:,1,batch_index)';
q(1,:) = bdb.data.q(:,1,batch_index)';
% note u(1,:) is the first control input and must be calculated

%%%% Handle input constraints and phases %%%%

%%%% WE WON"T NEED THIS YOU CAN TAKE IT OUT %%%%
% tphase = [0 0.5 2.5];  % the phase transition times [hours] (First element should be start of batch)
% Input constraints (columns for phase, rows for varia bles):

umin = [290];
umax = [323];

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

% weight_0 = .000005;
% weight_f = .000005;
% trans_start = 150;
% trans_end = 195;
% weight = ;

% for i = 1:bdb.nFE+1
%    if(i<trans_start)
%        input_move_weights(i) = weight_0;
%    elseif (i<trans_end)
%        input_move_weights(i) = weight_0 + ( weight_f-weight_0)*(cos(-pi + pi*(i-trans_start)/(trans_end - trans_start))+1)/2;
%    else
%        input_move_weights(i) = weight_f;
%    end 
% end

input_move_weights = weight;


for i = 1:bdb.nFE
    num_initial = size(Q_model.x0_index,2);
    num_meas = size(Q_model.x_index,2)+size(Q_model.y_index,2) + size(Q_model.u_index,2);
    %u(i:end,:) = bdb.unom(:,i:end)'%%% TEMPORARY REMOVE THIS CODE %%%
    
    %%%% create gams sets %%%%
    % time steps to batch termination T
    T.name = 'T';
    for j = 1:bdb.nFE+2-i
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
    data = data - repmat(Mult_model.Measure_Means, bdb.nFE+1-(i-1),1);
    data = data./ repmat(Mult_model.Measure_Std, bdb.nFE+1- (i-1),1);
    
    
    data = wgdxvar('data', data, 2);

    %%%% Determine Input Constraints %%%%
    % get intput means and stds
    input_meas_means = Mult_model.Measure_Means(size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)...
        +1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2) ...
        + size(Mult_model.Measurements{u_typ_idx},2));
    input_meas_std = Mult_model.Measure_Std(size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)...
        +1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2) ...
        + size(Mult_model.Measurements{u_typ_idx},2));
    
    for j = i:bdb.nFE+1
        phase = find(tphase<=j*bdb.samt, 1,'last');
        input_const(:,:,j-i+1) = [umin(:,phase) umax(:,phase)];
        
        % scale constraint
        input_const(:,:,j-i+1) = (input_const(:,:,j-i+1) - repmat(input_meas_means',1,2))./repmat(input_meas_std',1,2);
    end
    

    input_const = wgdxvar('input_const', input_const, 3);
    
    
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

    for j = i:bdb.nFE+1
       qual_coeffs(:,:,j-i+1) = Q_model.coeffs(num_initial+(j-1)*num_meas+1:num_initial+j*num_meas,:);
    end
    
    qaul_coeffs = wgdxvar('qual_coeffs', qual_coeffs, 3);
    
    %%%% Get Quality Model Scaling %%%%
    for j = i:bdb.nFE+1
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
    
    for j = i:bdb.nFE+1
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
    
    input_move_weights_wgdx = wgdxvar('input_move_weights', input_move_weights(1,i:end), 1);
    
    
    wgdx('MPCinitialize', data, mult_coeffs, mult_coeffs_intercept,PLS_intercept, ...
        Clust_Cent, meas_mean, meas_std, pred_mean, pred_std, contrib_past,...
        qaul_coeffs, q_meas_mean, q_meas_std, q_pred_mean, q_pred_std, qual_set, input_const, std_scores,...
        t_past_comp_gdx, w_star_reformed, Tsqrd_limit, past_data_gdx, w_star_past, SPE_limit, qual_lo, input_move_weights_wgdx);
    
    
    
    
    tic
    [x_out JMPC_out q_pred , ~, Tsqrd_out SPE_out] = gams('QBMPC_PMMA_gams nlp=conopt lo=0');
    time(i) = toc;
    
    if(mod(i,10) == 0)
        fprintf(['*** Batch ' num2str(batch_index) ' ***\n'])
    end
    fprintf(['  i = ' num2str(i) ' JMPC = ' num2str(JMPC_out.val) '\n']);
    
    JMPC(i) = JMPC_out.val;
    
    x_out = x_out.val;
    q_pred = q_pred.val;
    x_out = x_out.*repmat(Mult_model.Measure_Std,bdb.nFE-i+2,1) + repmat(Mult_model.Measure_Means,bdb.nFE-i+2,1);
    q_pred = q_pred'.*Q_model.y_std + Q_model.y_mean;
    
    q_pred_traj(i,:) = q_pred(:);
    
    % extract input moves (for the rest of the batch) from the gams output
    u(i:end-1,:) = x_out(1:end-1,size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2)...
        +1:size(Mult_model.Measurements{x_typ_idx},2)+size(Mult_model.Measurements{y_typ_idx},2) ...
        + size(Mult_model.Measurements{u_typ_idx},2));
    
    u_all(i:end,:,i) = u(i:end-1,:);
    
    
    %%%% PLOT INPUT TRAJECTORIES %%%%
  %  subplot(2,1,1)
    plot(i:241, x_out(:,4));
    hold on;
    plot(batch_data.unom(1,:),'r');
    plot(u(1:i,1), 'g', 'LineWidth', 1.5);
    
    %check SPE and T-squared
    Qstds = Q_model.x_std(num_initial+(i-1)*num_meas+1:end);
    Qmeans =Q_model.x_mean(num_initial+(i-1)*num_meas+1:end);
    batch_trajs_unfolded = [past_data...
        (reshape(x_out',1,size(x_out,1)*size(x_out,2))-Qmeans)./Qstds];
    scores = batch_trajs_unfolded*Q_model.wstar;
    
    Tsqrd(i) = sum((scores./Q_model.sa).^2);
    
    x_hat = scores*Q_model.wstar';
    error = batch_trajs_unfolded-x_hat;
    
    SPE(i) = error*error';
    
    %%%% Simulate to end of batch with the given input trajectory
    predict_to_end = false;
    if(predict_to_end)
        fprintf(['--- Simulating batch to completion \n']);
        X_out_complete(1,:) = x(i,:);
        T_out_complete = 0;
        for j = i:bdb.nFE
            [T_out X_out] = ode23s(@methyl_methacrylate, [(j-1)*bdb.samt j*bdb.samt], X_out_complete(end,:), '', x(1,:), u_all(j,:,i) , 'state');
            T_out_complete = [T_out_complete; T_out;];
            X_out_complete = [X_out_complete; X_out;];
            aux = methyl_methacrylate((i)*bdb.samt,X_out_complete(end,:), x(1,:), u_all(j,:,i), 'output');
            y_act(j+1,:,i) = [X_out_complete(end,9) aux(5:6) log(aux(5)) aux(7)]; %+ noise(i+1,:);
        end
        aux = methyl_methacrylate((i)*bdb.samt,X_out_complete(end,:), x(1,:), u_all(end,:,i), 'output');
        q_act(i,:) = aux(1:4);
    end
    
    [T_out X_out] = ode23s(@methyl_methacrylate, [(i-1)*bdb.samt i*bdb.samt], x(i,:), '', x(1,:), u(i,1) , 'state');
    x(i+1,:) = X_out(end,:);

    aux = methyl_methacrylate((i)*bdb.samt,x(i+1,:), x(1,:), u(i,1), 'output');
    y(i+1,:) = [x(i+1,9) aux(5:6) log(aux(5)) aux(7)]; %+ noise(i+1,:);
    q(i+1,:) = aux(1:4);
    
    %%%% Epty GAMs variables %%%%
    clear data contrib_past qaul_coeffs q_meas_mean q_meas_std q_pred_mean q_pred_std qual_set input_const w_star_reformed std_scores Tsqrd_limit t_past_comp;
    
    %%%% Empty GAMs Sets %%%%
    clear T M P C Q U max_min;
    

end
quality_result = q(end,1:3);

q_old = batch_data.q(1:3,end,batch_index)';

percent_imp = (quality_result - q_old)./(qual_target - q_old) * 100;
for i = 1:size(percent_imp,2)
   if(percent_imp(i)>100)
       percent_imp(i) = 200 - percent_imp(i);
   end
end
handle = gcf;

figure;
subplot(3,1,1)
plot(q_pred_traj(:,1));
hold on
plot([0 size(q_pred_traj,1)], [qual_target(1) qual_target(1)], '--k')
if(predict_to_end)
   plot(q_act(:,1),'r'); 
end

subplot(3,1,2)
plot(q_pred_traj(:,2));
hold on
plot([0 size(q_pred_traj,1)], [qual_target(2) qual_target(2)], '--k')
if(predict_to_end)
   plot(q_act(:,2),'r'); 
end

subplot(3,1,3)
plot(q_pred_traj(:,3));
hold on
plot([0 size(q_pred_traj,1)], [qual_target(3) qual_target(3)], '--k')
if(predict_to_end)
   plot(q_act(:,3),'r'); 
end

figure;
subplot(2,1,1)
plot(Tsqrd)
hold on;
plot([0 bdb.nFE+1], [Q_model.T95 Q_model.T95], ':k')
plot([0 bdb.nFE+1], [Q_model.T99 Q_model.T99], '--k')
xlabel('Sampling instant');
ylabel('T^2');

subplot(2,1,2)
plot(SPE)
hold on;
plot([0 bdb.nFE+1], [Q_model.SPE95 Q_model.SPE95], ':k')
plot([0 bdb.nFE+1], [Q_model.SPE99 Q_model.SPE99], '--k')
xlabel('Sampling instant');
ylabel('SPE');




end