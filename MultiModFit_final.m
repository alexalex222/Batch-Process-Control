clear variables;
close all;

load('test4.mat');

% This is the dynmaic model for sample 401-600

IndexQinter1=201;
IndexQinter2=401;
PredictionStep=200;

nb_fit = 45; % number of batches to use fitting MultiMod_finalel
nb_val = 15; % number of batches to use validating MultiMod_finalel



u_fit = zeros(size(bdb.data.u,1), PredictionStep+1, nb_fit);
u_val = zeros(size(bdb.data.u,1), PredictionStep+1, nb_val);

u_fit(:,1:end-1,1:nb_fit-2) = bdb.data.u(:,IndexQinter2:end,1:nb_fit-2);
u_fit(:,1:end-1,nb_fit-1:nb_fit) = bdb.data.u_PRBS(:,IndexQinter2:end,1:2);

u_val(:,1:end-1,:) = bdb.data.u(:,IndexQinter2:end,nb_fit+1:nb_fit+nb_val);

y_fit = zeros(size(bdb.data.y,1), PredictionStep+1, nb_fit);
y_val = zeros(size(bdb.data.y,1), PredictionStep+1, nb_val);


% I add two batches of PRBS into y_fit
y_fit(:,:,1:nb_fit-2) = bdb.data.y(:,IndexQinter2:end,1:nb_fit-2);
y_fit(:,:,nb_fit-1:nb_fit) = bdb.data.y_PRBS(:,IndexQinter2:end,1:2);
y_val = bdb.data.y(:,IndexQinter2:end,nb_fit+1:nb_fit+nb_val);


% "aligning inputs" - so that the input vector is as long as the state
for i = 1:size(u_fit,3)
   u_fit(:,end,i) =  u_fit(:,end-1,i); 
end

for i = 1:size(u_val,3)
   u_val(:,end,i) = u_val(:,end-1,i);
end

% u_example=u_fit(:,:,1);
% u_example=reshape(u_example,6,601);

% Fit a mult MultiMod_finalel to the data
MultiMod_final = Mult_Mod();
MultiMod_final.SetRegressionType('OLS'); % This can be 'PLS or OLS' if it is PLS you need to "play" with the PLSnPCs
MultiMod_final.SetPLSnPC(10); % THIS CAN BE A VECTOR BUT FOR NOW IT WON"T MATTER
MultiMod_final.SetClusters(2:5); % THIS is a key parameter to vary
%MultiMod_final.SetPLSnPC(24);
%MultiMod_final.SetClusters(5);
MultiMod_final.SetBatchNumber(nb_fit, nb_val);
MultiMod_final.SetBatchLength(PredictionStep+1, bdb.samt);
%    MultiMod_final.Cluster_Iter = 15;

% THiS SETS UP THE MultiMod_finalEL STRUCTURE

% MultiMod_final.addMeasureed('type', index  in the database, number of lags)
MultiMod_final.AddMeasured('input', 1,1);
MultiMod_final.AddMeasured('input', 2,1);
MultiMod_final.AddMeasured('input', 3,1);
% MultiMod_final.AddMeasured('state',1,1); % note that we'll be using outputs for this system not states
% MultiMod_final.AddMeasured('state',2,1);
% MultiMod_final.AddMeasured('state',3,1);

MultiMod_final.AddMeasured('output',2,1); % note that we'll be using outputs for this system not states
MultiMod_final.AddMeasured('output',3,1);
MultiMod_final.AddMeasured('output',4,1);
MultiMod_final.AddMeasured('output',5,1);

% MultiMod_final.AddPredicted('state',1);
% MultiMod_final.AddPredicted('state',2);
% MultiMod_final.AddPredicted('state',3);

MultiMod_final.AddPredicted('output',2);
MultiMod_final.AddPredicted('output',3);
MultiMod_final.AddPredicted('output',4);
MultiMod_final.AddPredicted('output',5);



MultiMod_final.Fit_Model('input_fit', u_fit,'input_val', u_val,...
    'output_fit', y_fit,'output_val',y_val); % this should be output_fit rather then state_fit

save('MultiMod_final','MultiMod_final');

for i=1:5
    MultiMod_final.Plot_Traj('input', u_val(:,:,i),...
        'output',y_val(:,:,i));
end