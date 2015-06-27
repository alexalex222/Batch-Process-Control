clear variables;
close all;

load('test4.mat');

nb_fit = 45; % number of batches to use fitting model
nb_val = 15; % number of batches to use validating model

% get data out of the database
% x_fit = bdb.data.x(:,:,1:nb_fit);
% x_val = bdb.data.x(:,:,nb_fit+1:nb_fit+nb_val);

u_fit = zeros(size(bdb.data.u,1), size(bdb.data.u,2)+1, nb_fit);
u_val = zeros(size(bdb.data.u,1), size(bdb.data.u,2)+1, nb_val);

u_fit(:,1:end-1,1:nb_fit-2) = bdb.data.u(:,:,1:nb_fit-2);
u_fit(:,1:end-1,nb_fit-1:nb_fit) = bdb.data.u_PRBS(:,:,1:2);

u_val(:,1:end-1,:) = bdb.data.u(:,:,nb_fit+1:nb_fit+nb_val);

y_fit = zeros(size(bdb.data.y,1), size(bdb.data.y,2), nb_fit);
y_val = zeros(size(bdb.data.y,1), size(bdb.data.y,2), nb_val);


% I add two batches of PRBS into y_fit
y_fit(:,:,1:nb_fit-2) = bdb.data.y(:,:,1:nb_fit-2);
y_fit(:,:,nb_fit-1:nb_fit) = bdb.data.y_PRBS(:,:,1:2);
y_val = bdb.data.y(:,:,nb_fit+1:nb_fit+nb_val);


% "aligning inputs" - so that the input vector is as long as the state
for i = 1:size(u_fit,3)
   u_fit(:,end,i) =  u_fit(:,end-1,i); 
end

for i = 1:size(u_val,3)
   u_val(:,end,i) = u_val(:,end-1,i);
end

% u_example=u_fit(:,:,1);
% u_example=reshape(u_example,6,601);

% Fit a mult model to the data
mod = Mult_Mod();
mod.SetRegressionType('OLS'); % This can be 'PLS or OLS' if it is PLS you need to "play" with the PLSnPCs
mod.SetPLSnPC(10); % THIS CAN BE A VECTOR BUT FOR NOW IT WON"T MATTER
mod.SetClusters(2); % THIS is a key parameter to vary
%mod.SetPLSnPC(24);
%mod.SetClusters(5);
mod.SetBatchNumber(nb_fit, nb_val);
mod.SetBatchLength(bdb.nFE+1, bdb.samt);
%    mod.Cluster_Iter = 15;

% THiS SETS UP THE MODEL STRUCTURE

% mod.addMeasureed('type', index  in the database, number of lags)
mod.AddMeasured('input', 1,1);
mod.AddMeasured('input', 2,1);
mod.AddMeasured('input', 3,1);
mod.AddMeasured('input', 4,1);

% mod.AddMeasured('state',1,1); % note that we'll be using outputs for this system not states
% mod.AddMeasured('state',2,1);
% mod.AddMeasured('state',3,1);

mod.AddMeasured('output',2,1); % note that we'll be using outputs for this system not states
mod.AddMeasured('output',3,1);
mod.AddMeasured('output',4,1);
mod.AddMeasured('output',5,1);

% mod.AddPredicted('state',1);
% mod.AddPredicted('state',2);
% mod.AddPredicted('state',3);

mod.AddPredicted('output',2);
mod.AddPredicted('output',3);
mod.AddPredicted('output',4);
mod.AddPredicted('output',5);

% mod.Fit_Model('input_fit', u_fit,'input_val', u_val,...
%     'state_fit', x_fit,'state_val',x_val); % this should be output_fit rather then state_fit

mod.Fit_Model('input_fit', u_fit,'input_val', u_val,...
    'output_fit', y_fit,'output_val',y_val); % this should be output_fit rather then state_fit

save('mod1','mod');

for i=1:5
    mod.Plot_Traj('input', u_val(:,:,i),...
        'output',y_val(:,:,i));
end