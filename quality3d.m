clear;
clc;
close all;


load PenMPCresults;

qual_target = [0.0155 12.9297 1.3280];
% qual_target =mean(Quality_MPC);

figure;
scatter3(Quality_PI(1:20,1),Quality_PI(1:20,2),Quality_PI(1:20,3),'.');
hold on;
scatter3(Quality_MPC(1:20,1),Quality_MPC(1:20,2),Quality_MPC(1:20,3));

scatter3(qual_target(1,1),qual_target(1,2),qual_target(1,3),100,'rp');
legend('PI controller', 'QBMPC', 'Desired Quality');
xlabel('Glucose g/L');
ylabel('Biomass g/L');
zlabel('Peniclllin g/L');