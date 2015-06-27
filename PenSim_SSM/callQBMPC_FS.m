clear;
clc;
close all;

load SSM_Quality;


load test9;

Num_Batch=60; % Number of batches to be tested in the MPC



Quality_MPC=zeros(Num_Batch, bdb.nq);
Quality_PI=zeros(Num_Batch, bdb.nq);

for i=1:Num_Batch
    [quality_result, q_old, x, y,q, u, JMPC, time] =...
          QBMPC_FinalState(bdb, i, SSM);
    Quality_MPC(i,:) = quality_result;
    Quality_PI(i,:)= q_old;
end

save('PenMPCresults.mat','Quality_MPC','Quality_PI');