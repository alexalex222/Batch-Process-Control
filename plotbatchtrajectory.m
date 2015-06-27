clear;
clc;
close all;

load u_MPC;
load y_MPC;
load u_MPC_inter;
load y_MPC_inter;
load test9;
load Y1;
load Y2;
load Y3;


offset1= [0: 4.1667e-04:0.25];
offset1=offset1(1:600)';

u_PI = bdb.data.u(:,:,7);
u_PI = u_PI';
u_PI(:,4) = u_PI(:,4)+offset1;

y_PI = bdb.data.y(2:5,:,7);
y_PI = y_PI';

y_MPC = y_MPC(:,2:5);
y_MPC_inter = y_MPC_inter(:,2:5);

figure;
subplot 221;
stairs(u_PI(:,1));
hold on;
stairs(u_MPC(:,1),'k');
stairs(u_MPC_inter(:,1),'r');
xlabel('Sample');
ylabel('Feed Rate (L/h)');
hold off

subplot 222;
stairs(u_PI(:,2));
hold on;
stairs(u_MPC(:,2),'k');
stairs(u_MPC_inter(:,2),'r');
xlabel('Sample');
ylabel('Aeration Rate (L/h)');
hold off
legend('PI','Single Model based MPC', 'Multiple Models based MPC');

subplot 223;
stairs(u_PI(:,3));
hold on;
stairs(u_MPC(:,3),'k');
stairs(u_MPC_inter(:,3),'r');
xlabel('Sample');
ylabel('Agitator Power (w)');
hold off

subplot 224;
stairs(u_PI(:,4));
hold on;
stairs(u_MPC(:,4),'k');
stairs(u_MPC_inter(:,4),'r');
xlabel('Sample');
ylabel('Cooling Water Flow Rate (L/H)');
hold off

figure;
subplot 221;
stairs(Y1(:,2));
hold on;
stairs(Y2(:,2),'k');
stairs(Y3(:,2),'r');
xlabel('Sample');
ylabel('Reactor Temperature (K)');
xlim([0 600]);
hold off

subplot 222;
stairs(y_PI(:,2));
hold on;
stairs(y_MPC(:,2),'k');
stairs(y_MPC_inter(:,2),'r');
xlabel('Sample');
ylabel('Dissoved O_2 (g/L)');
xlim([0 600]);
hold off
legend('PI','Single Model based MPC', 'Multiple Models based MPC');

subplot 223;
stairs(y_PI(:,3));
hold on;
stairs(y_MPC(:,3),'k');
stairs(y_MPC_inter(:,3),'r');
xlabel('Sample');
ylabel('Dissoved CO_2 (g/L)');
xlim([0 600]);
hold off

subplot 224;
stairs(y_PI(:,4));
hold on;
stairs(y_MPC(:,4),'k');
stairs(y_MPC_inter(:,4),'r');
xlabel('Sample');
ylabel('Fermentor Volumn (L)');
xlim([0 600]);
hold off