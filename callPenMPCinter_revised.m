clear;
clc;
close all;

load Q_model_Inter1;
load Q_model_Inter2
load Q_model_final;
load mod1;


load test4;

Num_Batch=50; % Number of batches to be tested in the MPC
Quality_MPC_inter=zeros(Num_Batch, bdb.nq);


for i=1:Num_Batch
    [quality_result, q_old, x, y,q, u, JMPC, time] =...
          QBMPC_Pen_Inter_Re(bdb, i, Q_model_Inter1, Q_model_Inter2, Q_model_Final,...
          mod);
    Quality_MPC_inter(i,:) = quality_result;
    
end

save('PenMPCinter_Re2','Quality_MPC_inter');

% load Q_model_v1;
% load mod1;



% Quality_MPC=zeros(Num_Batch, bdb.nq);
% Quality_PI=zeros(Num_Batch, bdb.nq);
% 
% for i=1:Num_Batch
%     [quality_result, q_old, x, y,q, u, JMPC, time] =...
%           QBMPC_PenV1(bdb, i, Q_model, mod);
%     Quality_MPC(i,:) = quality_result;
%     Quality_PI(i,:)= q_old;
% end
% 
% save('PenMPCresults3.mat','Quality_MPC');