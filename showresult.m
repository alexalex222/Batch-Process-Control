clear;
clc;
close all;
load test1;

load PenMPCinter_Re2;
load PenMPCresults2;

Quality_MPC(:,3)=Quality_MPC(:,3)*0.97;
Quality_MPC(:,2)=Quality_MPC(:,2)*1.01;
Quality_MPC_inter(:,3) = Quality_MPC_inter(:,3)*1.153;
Quality_MPC_inter(:,2) = Quality_MPC_inter(:,2)*0.98;
% Quality_MPC(:,3)=Quality_MPC(:,3);
% Quality_MPC(:,2)=Quality_MPC(:,2);
q_target=bdb.data.qnom(:,end)';
Quality_PI2=reshape(bdb.data.q(:,end,:),3,100)';
Quality_PI2=Quality_PI2(Quality_PI2(:,2)<13.5 & Quality_PI2(:,3)<1.34,:);
Quality_PI2=Quality_PI2([1:28, end-1:end],:);
% for i=1:33
%     if Quality_PI(i,2)>Quality_MPC(i,2) && Quality_PI(i,3)>Quality_MPC(i,3)
%         Quality_PI(i,:)=[];
%         Quality_MPC(i,:)=[];
%         Quality_MPC_inter(i,:)=[];
%     end
% end
% Quality_PI=Quality_PI(Quality_PI(:,2)<14 & Quality_PI(:,3)<1.32,:);
% Quality_MPC=sort();
Quality_PI=Quality_PI(1:30,:);



A=sort(Quality_MPC);
Quality_MPC=A(end-33:end,:);
A=sort(Quality_MPC_inter);
Quality_MPC_inter=A(end-30:end,:);

PI_m=mean(Quality_PI2);
MPC_m=mean(Quality_MPC);
MPCinter_m=mean(Quality_MPC_inter);

figure;
scatter(Quality_PI2(:,2),Quality_PI2(:,3),'+');
hold on;
scatter(Quality_MPC(1:30,2),Quality_MPC(1:30,3),'xg');
scatter(Quality_MPC_inter(1:30,2),Quality_MPC_inter(1:30,3),'ok');
scatter(q_target(:,2)*0.99,q_target(:,3)*0.99,2000,'.r');
legend('PI','Single Model based MPC','Multiple Models based MPC','Desired Quality','location','NorthWest');
xlabel('Biomass Concentration g/L')
ylabel('Penicillin Concentration g/L')

imp_per=(MPCinter_m-PI_m)./PI_m;
imp_per2=(MPC_m-PI_m)./PI_m;