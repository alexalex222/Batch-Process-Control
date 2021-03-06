clear;
clc;
close all;

load test1;
load SSM_Quality;

A = SSM.A;
B = SSM.B;
C = SSM.C;
D = SSM.D;
J_k = SSM.J_k;
q_meas_mean = SSM.mt_final_mean;
q_meas_std = SSM.mt_final_std;
q_pred_mean = SSM.q_final_mean;
q_pred_std = SSM.q_final_std;
FuureHorizon = SSM.FutureHorizon;
PastHorizon = SSM.FutureHorizon;
% k = SSM.k;
Coeff_Quality = SSM.Coeff_Quality;
y_typ_idx = SSM.y_typ_idx;
u_typ_idx = SSM.u_typ_idx;

% select a test batch
for i = 19: 19;

    Y = bdb.data.y(y_typ_idx,:,i);
    U = bdb.data.u(u_typ_idx,:,i);
    U(:,end+1) = U(:,end);

    mt_initial = [bdb.data.y(y_typ_idx,1,i)'...
        bdb.data.y(y_typ_idx,2,i)' bdb.data.u(u_typ_idx,1,i)' bdb.data.u(u_typ_idx,2,i)']*J_k';
    
    
    Y_hat =zeros(size(Y));

    % set the current sampling instance
    k = 3;
    Y_hat(:,1:k-1) = Y(:,1:k-1);

    for j = k:601
        pt = [Y_hat(:,j-2)' Y_hat(:,j-1)' U(:,j-2)'+[normrnd(0,0.02), normrnd(0,1), normrnd(0,2), normrnd(0,0.01) ] U(:,j-1)'];
%       pt = [ Y_hat(:,j-1)'  U(:,j-1)'];
        Y_hat(:,j) = SSM.C*SSM.J_k*pt';
    end  



% mt_final_hat = [Y_hat(:,end-2)' Y_hat(:,end-1)' U(:,end-2)' U(:,end-1)']*J_k';

mt_final_hat = [Y_hat(:,end-1)' Y_hat(:,end)' U(:,end-1)' U(:,end)']*J_k';

if SSM.ScaldQualityModel ==1
    q_final_hat = (mt_final_hat-q_meas_mean)./q_meas_std*Coeff_Quality(2:end,:).*q_pred_std+q_pred_mean;
else
    q_final_hat = mt_final_hat*Coeff_Quality(1:end,:);
end




%     data = [Y_hat(:,k:600)' U(:,k:600)'];
% 
%     Coeff_Quality_temp = Convert3DCoeff_Quality(Coeff_Quality, J_k, PastHorizon, 6, 4);
% 
%     [d1,~] = size(data);
%     Coeff_Quality = zeros(size(SSM.y_typ_idx,2)+ size(SSM.u_typ_idx,2),bdb.nq,d1);
%     Coeff_Quality(:,:,end-1:end) = Coeff_Quality_temp;
% 
%     q_final_hat = zeros(3,1);   
%     
%     for Q = 1:3
%         for T = 1:d1
%             for M = 1:size(data,2)
%             
%                 q_final_hat(Q) = q_final_hat(Q)+data(T,M).*Coeff_Quality(M,Q,T);
%             
%             end
%         end 
%     
%     end

    q_final_hat_traj(i,:) = q_final_hat';

    q_final_traj(i,:) = bdb.data.q(:,end,i)';



end

%% plot the reuslts

Y = Y';
Y_hat = Y_hat';

figure;
subplot 321
plot(Y(:,1));
hold on
plot(Y_hat(:,1), 'r--');
xlim([0 600]);
title('Jacket Temperature');
xlabel('Time');
ylabel('K');
legend('Actual','Estimated');
R2y1=rsquare(Y_hat(:,1),Y(:,1));

subplot 322
plot(Y(:,2));
hold on
plot(Y_hat(:,2), 'r--');
xlim([0 600]);
title('Fermentor Temperature');
xlabel('Time');
ylabel('K');
R2y2=rsquare(Y_hat(:,2),Y(:,2));

subplot 323
plot(Y(:,3));
hold on
plot(Y_hat(:,3), 'r--');
xlim([0 600]);
title('Dissolved O_2');
xlabel('Time');
ylabel('mmol/L');
R2y3=rsquare(Y_hat(:,3),Y(:,3));

subplot 324
plot(Y(:,4));
hold on
plot(Y_hat(:,4), 'r--');
xlim([0 600]);
title('CO_2');
xlabel('Time');
ylabel('mmol/L');
R2y4=rsquare(Y_hat(:,4),Y(:,4));

subplot 325
plot(Y(:,5));
hold on
plot(Y_hat(:,5), 'r--');
xlim([0 600]);
title('Fermentor Volume');
xlabel('Time');
ylabel('L');
R2y5=rsquare(Y_hat(:,5),Y(:,5));

subplot 326
plot(Y(:,6));
hold on
plot(Y_hat(:,6), 'r--');
xlim([0 600]);
title('Feed Temperature');
xlabel('Time');
ylabel('K');
R2y6=rsquare(Y_hat(:,6),Y(:,6));


figure;
subplot 311
plot(q_final_traj(:,1),'o');
hold on;
plot(q_final_hat_traj(:,1),'xr');
title('q_1');
legend('Actual','Predicted');


subplot 312
plot(q_final_traj(:,2),'o');
hold on;
plot(q_final_hat_traj(:,2),'xr');
title('q_1');


subplot 313
plot(q_final_traj(:,3),'o');
hold on;
plot(q_final_hat_traj(:,3),'xr');
title('q_3');


