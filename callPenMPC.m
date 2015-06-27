clear;
clc;
close all;

load Q_model_v1;
load mod;
load test9;

Num_Batch=50; % Number of batches to be tested in the MPC
Quality_MPC=zeros(Num_Batch, bdb.nq);
Quality_PI=zeros(Num_Batch, bdb.nq);

for i=1:Num_Batch
    [quality_result, q_old, x, y,q, u, JMPC, time] =...
          QBMPC_PenV1(bdb, i, Q_model, mod);
    Quality_MPC(i,:) = quality_result;
    Quality_PI(i,:)= q_old;
end

save('PenMPCresults.mat','Quality_MPC','Quality_PI');