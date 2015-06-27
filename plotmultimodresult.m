clear variables;
close all;

load('test4.mat');

% This is the dynmaic model for sample 401-600

IndexQinter1=201;
IndexQinter2=401;
PredictionStep=200;

nb_fit = 45; % number of batches to use fitting MultiMod_inter1el
nb_val = 15; % number of batches to use validating MultiMod_inter1el



u_fit = zeros(size(bdb.data.u,1), PredictionStep+1, nb_fit);
u_val = zeros(size(bdb.data.u,1), PredictionStep+1, nb_val);

u_fit(:,1:end-1,1:nb_fit) = bdb.data.u(:,1:IndexQinter1-1,1:nb_fit);
% u_fit(:,1:end-1,nb_fit-1:nb_fit) = bdb.data.u_PRBS(:,1:IndexQinter1-1,1:2);

u_val(:,1:end-1,:) = bdb.data.u(:,1:IndexQinter1-1,nb_fit+1:nb_fit+nb_val);

y_fit = zeros(size(bdb.data.y,1), PredictionStep+1, nb_fit);
y_val = zeros(size(bdb.data.y,1), PredictionStep+1, nb_val);


% I add two batches of PRBS into y_fit
y_fit(:,:,1:nb_fit) = bdb.data.y(:,1:IndexQinter1,1:nb_fit);
% y_fit(:,:,nb_fit-1:nb_fit) = bdb.data.y_PRBS(:,1:IndexQinter1,1:2);
y_val = bdb.data.y(:,1:IndexQinter1,nb_fit+1:nb_fit+nb_val);


% "aligning inputs" - so that the input vector is as long as the state
for i = 1:size(u_fit,3)
   u_fit(:,end,i) =  u_fit(:,end-1,i); 
end

for i = 1:size(u_val,3)
   u_val(:,end,i) = u_val(:,end-1,i);
end

% u_example=u_fit(:,:,1);
% u_example=reshape(u_example,6,601);

load MultiMod_inter1;

figure

A=MultiMod_inter1.Plot_Traj('input', u_val(:,:,9),...
        'output',y_val(:,:,9));



% A=MultiMod_inter1.Plot_Traj('input', u_val(:,:,6),'output',y_val(:,:,6));


clear;

load('test4.mat');

% This is the dynmaic model for sample 401-600

IndexQinter1=201;
IndexQinter2=401;
PredictionStep=200;

nb_fit = 45; % number of batches to use fitting MultiMod_inter2el
nb_val = 15; % number of batches to use validating MultiMod_inter2el



u_fit = zeros(size(bdb.data.u,1), PredictionStep+1, nb_fit);
u_val = zeros(size(bdb.data.u,1), PredictionStep+1, nb_val);

u_fit(:,1:end-1,1:nb_fit) = bdb.data.u(:,IndexQinter1:IndexQinter2-1,1:nb_fit);
% u_fit(:,1:end-1,nb_fit-1:nb_fit) = bdb.data.u_PRBS(:,IndexQinter1:IndexQinter2-1,1:2);

u_val(:,1:end-1,:) = bdb.data.u(:,IndexQinter1:IndexQinter2-1,nb_fit+1:nb_fit+nb_val);

y_fit = zeros(size(bdb.data.y,1), PredictionStep+1, nb_fit);
y_val = zeros(size(bdb.data.y,1), PredictionStep+1, nb_val);


% I add two batches of PRBS into y_fit
y_fit(:,:,1:nb_fit) = bdb.data.y(:,IndexQinter1:IndexQinter2,1:nb_fit);
% y_fit(:,:,nb_fit-1:nb_fit) = bdb.data.y_PRBS(:,IndexQinter1:IndexQinter2,1:2);
y_val = bdb.data.y(:,IndexQinter1:IndexQinter2,nb_fit+1:nb_fit+nb_val);


% "aligning inputs" - so that the input vector is as long as the state
for i = 1:size(u_fit,3)
   u_fit(:,end,i) =  u_fit(:,end-1,i); 
end

for i = 1:size(u_val,3)
   u_val(:,end,i) = u_val(:,end-1,i);
end
load MultiMod_inter2;

figure
% for i=[1 2 3 4 5]
A=MultiMod_inter2.Plot_Traj('input', u_val(:,:,9),...
        'output',y_val(:,:,9));
% end

clear;

load('test4.mat');

% This is the dynmaic model for sample 401-600

IndexQinter1=201;
IndexQinter2=401;
PredictionStep=200;

nb_fit = 45; % number of batches to use fitting MultiMod_finalel
nb_val = 15; % number of batches to use validating MultiMod_finalel



u_fit = zeros(size(bdb.data.u,1), PredictionStep+1, nb_fit);
u_val = zeros(size(bdb.data.u,1), PredictionStep+1, nb_val);

u_fit(:,1:end-1,1:nb_fit) = bdb.data.u(:,IndexQinter2:end,1:nb_fit);
% u_fit(:,1:end-1,nb_fit-1:nb_fit) = bdb.data.u_PRBS(:,IndexQinter2:end,1:2);

u_val(:,1:end-1,:) = bdb.data.u(:,IndexQinter2:end,nb_fit+1:nb_fit+nb_val);

y_fit = zeros(size(bdb.data.y,1), PredictionStep+1, nb_fit);
y_val = zeros(size(bdb.data.y,1), PredictionStep+1, nb_val);


% I add two batches of PRBS into y_fit
y_fit(:,:,1:nb_fit) = bdb.data.y(:,IndexQinter2:end,1:nb_fit);
% y_fit(:,:,nb_fit-1:nb_fit) = bdb.data.y_PRBS(:,IndexQinter2:end,1:2);
y_val = bdb.data.y(:,IndexQinter2:end,nb_fit+1:nb_fit+nb_val);


% "aligning inputs" - so that the input vector is as long as the state
for i = 1:size(u_fit,3)
   u_fit(:,end,i) =  u_fit(:,end-1,i); 
end

for i = 1:size(u_val,3)
   u_val(:,end,i) = u_val(:,end-1,i);
end

load MultiMod_final;

figure
% for i=[1 2 3 4 5]
    A=MultiMod_final.Plot_Traj('input', u_val(:,:,4),...
        'output',y_val(:,:,4));
% end