clear;
clc;
close all;
load PenMPCresults;
load test9;

load U_PI;
load U_MPC;

q_target=bdb.data.qnom(:,end)';

Quality_MPC = Quality_MPC(1:30,:);
Quality_PI = Quality_PI(1:50,:);

Quality_PI=Quality_PI(Quality_PI(:,2)<14 & Quality_PI(:,3)<1.37,:);

Quality_PI=Quality_PI(1:30,:);



PI_m=mean(Quality_PI);
MPC_m=mean(Quality_MPC);

figure;
scatter(Quality_PI(:,2),Quality_PI(:,3),'+');
hold on;
scatter(Quality_MPC(:,2),Quality_MPC(:,3),'xr');
scatter(q_target(:,2)*0.99,q_target(:,3)*0.99,1000,'.k');
legend('PI','CVA based MPC','Desired Quality');
xlabel('Biomass Concentration g/L')
ylabel('Penicillin Concentration g/L')


figure
subplot 221
stairs(U_PI(:,1));
hold on;
stairs(U_MPC(:,1),'r');
xlabel('Glucose feed rate');
ylabel('L/h');
legend('PI','MPC');


subplot 222
stairs(U_PI(:,2));
hold on;
stairs(U_MPC(:,2),'r');
xlabel('Aeration rate');
ylabel('L/h');

subplot 223
stairs(U_PI(:,3));
hold on;
stairs(U_MPC(:,3),'r');
xlabel('Agitator power');
ylabel('V');

subplot 224
stairs(U_PI(:,4));
hold on;
stairs(U_MPC(:,4),'r');
xlabel('Cooling water flow rate');
ylabel('L/h');