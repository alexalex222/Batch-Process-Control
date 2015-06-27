clear;
clc;
close all;

load test2;

% 7, 20, 31, 35, 43 are bad batches


% x_PRBS = bdb.data.x_PRBS;
% y_PRBS = bdb.data.y_PRBS;
% u_PRBS = bdb.data.u_PRBS;

bdb.data.x_PRBS(:,:, [7 20 31 35 43]) = [];
bdb.data.y_PRBS(:,:, [7 20 31 35 43]) = [];
bdb.data.u_PRBS(:,:, [7 20 31 35 43]) = [];

save('test_PRBS.mat','bdb');