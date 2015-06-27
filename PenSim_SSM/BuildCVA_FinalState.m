clear
clc;
close all;


load test1;

PastHorizon = 2;
FutureHorizon = 2;

% number of test batches
NTesBatch = 30;

idx_inter1 = 201;
idx_inter2 = 401;

% dimension of the state
k = 11;

% Select the x, y and u to be used in the model
x_typ_idx = [1 2 3 4 5 6 7];
y_typ_idx = [1 2 3 4 5 6];
u_typ_idx = [1 2 3 4];
q_typ_idx = [1 2 3];

[ P, F, Q, Xu, Yu, Uu ] = UnfoldBatch( bdb.data.x(:,:,1:end-NTesBatch),...
    bdb.data.y(y_typ_idx,:,1:end-NTesBatch), bdb.data.u(u_typ_idx,:,1:end-NTesBatch), PastHorizon, FutureHorizon );

% [Pscld, Pmean, Psigma] = zscore(P);
% [Fscld, Fmean, Fsigma] = zscore(F);
% [Qscld, Qmean, Qsigma] = zscore(Q);
% [~, Ymean, Ysigma] = zscore(Yu);
% [~, Umean, Usigma] = zscore(Uu);


[ J, L, D, J_k, SigmaXX, SigmaYY, Omega ] = CVA( P, F, Q, k );

[ A_hat, B_hat, C_hat, D_hat ] = SysID( P, F, Q, k, Yu, Uu  );



% select the test batch
i = 10;

% unfold the test data
[ Pt, Ft, Qt, Xt, Yt, Ut ] = UnfoldBatch( bdb.data.x(:,:,i:i),...
    bdb.data.y(y_typ_idx,:,i:i), bdb.data.u(u_typ_idx,:,i:i), PastHorizon, FutureHorizon );


Y_hatID = Pt*J_k'*C_hat'; 



figure;
subplot 321
plot(Yt(1:end,1));
hold on
plot(Y_hatID(1:end,1), 'r--');
title('y_1');
legend('Actual Y','Estimated Y');
R2y1=rsquare(Y_hatID(:,1),Yt(:,1));

subplot 322
plot(Yt(1:end,2));
hold on
plot(Y_hatID(1:end,2), 'r--');
title('y_2');
R2y2=rsquare(Y_hatID(:,2),Yt(:,2));

subplot 323
plot(Yt(1:end,3));
hold on
plot(Y_hatID(1:end,3), 'r--');
title('y_3');
R2y3=rsquare(Y_hatID(:,3),Yt(:,3));

subplot 324
plot(Yt(1:end,4));
hold on
plot(Y_hatID(1:end,4), 'r--');
title('y_4');
R2y4=rsquare(Y_hatID(:,4),Yt(:,4));

subplot 325
plot(Yt(1:end,5));
hold on
plot(Y_hatID(1:end,5), 'r--');
title('y_5');
R2y5=rsquare(Y_hatID(:,5),Yt(:,5));

subplot 326
plot(Yt(1:end,6));
hold on
plot(Y_hatID(1:end,6), 'r--');
title('y_6');
R2y6=rsquare(Y_hatID(:,6),Yt(:,6));

% build the quality model

mt_final = zeros(bdb.nb-NTesBatch,k);

mt_inter1 = zeros(bdb.nb-NTesBatch,k);

mt_inter2 = zeros(bdb.nb-NTesBatch,k);

q_final = zeros(bdb.nb-NTesBatch,size(q_typ_idx,2));

q_inter1 = zeros(bdb.nb-NTesBatch,size(q_typ_idx,2));

q_inter2 = zeros(bdb.nb-NTesBatch,size(q_typ_idx,2));

for j = 1: bdb.nb-NTesBatch
    
    mt_final(j,:) = [bdb.data.y(y_typ_idx,end-2,j)'...
        bdb.data.y(y_typ_idx,end-1,j)' bdb.data.u(u_typ_idx,end-2,j)' bdb.data.u(u_typ_idx,end-1,j)']*J_k';
    
    mt_inter1(j,:) = [bdb.data.y(y_typ_idx,idx_inter1-2,j)'...
        bdb.data.y(y_typ_idx,idx_inter1-1,j)' bdb.data.u(u_typ_idx,idx_inter1-2,j)' bdb.data.u(u_typ_idx,idx_inter1-1,j)']*J_k';
    
     mt_inter2(j,:) = [bdb.data.y(y_typ_idx,idx_inter2-2,j)'...
        bdb.data.y(y_typ_idx,idx_inter2-1,j)' bdb.data.u(u_typ_idx,idx_inter2-2,j)' bdb.data.u(u_typ_idx,idx_inter2-1,j)']*J_k';

    
end

for j = 1: bdb.nb-NTesBatch
    
    q_final(j,:) = bdb.data.q(q_typ_idx,end,j)';
    
    q_inter1(j,:) = bdb.data.q(q_typ_idx,idx_inter1,j)';
    
    q_inter2(j,:) = bdb.data.q(q_typ_idx,idx_inter2,j)';

end

x_final = bdb.data.x(x_typ_idx,end,j)';

for j = 1: bdb.nb-NTesBatch
    
    x_final(j,:) = bdb.data.x(x_typ_idx,end,j)';

end

[mt_final_scld, mt_final_mean, mt_final_std] = zscore(mt_final);

[q_final_scld, q_final_mean, q_final_std] = zscore(q_final);

[x_final_scld, x_final_mean, x_final_std] = zscore(x_final);

ScaldQualityModel = false;

if ScaldQualityModel==1
    [~,~,~,~,Coeff_Quality,~,MSE] = plsregress(mt_final_scld, q_final_scld, size(mt_final,2), 'cv', bdb.nb-NTesBatch);
    q_final_hat = [ones(bdb.nb-NTesBatch, 1) mt_final_scld]*Coeff_Quality.*repmat(q_final_std, bdb.nb-NTesBatch, 1)...
        +repmat(q_final_mean, bdb.nb-NTesBatch, 1);
else
    Coeff_Quality = ols([mt_inter1; mt_inter2; mt_final], [q_inter1; q_inter2; q_final]);
    q_final_hat = mt_final*Coeff_Quality(1:end,:);
end

ScaldRotationModel = false;

if ScaldRotationModel == 1
    [~,~,~,~,Rotation,~,MSE] = plsregress(mt_final,x_final,size(mt_final,2), 'cv', bdb.nb-NTesBatch);
    x_final_hat = [ones(bdb.nb-NTesBatch, 1) mt_final_scld]*Rotation.*repmat(x_final_std, bdb.nb-NTesBatch, 1)...
        +repmat(x_final_mean, bdb.nb-NTesBatch, 1);
else
    Rotation = ols(mt_final, x_final);
    x_final_hat = mt_final*Rotation;
end



figure;
subplot 311
plot(q_final(:,1),'o');
hold on;
plot(q_final_hat(:,1),'xr');
title('q_1');
legend('Actual','Predicted');
R2q1=rsquare(q_final_hat(:,1),q_final(:,1));

subplot 312
plot(q_final(:,2),'o');
hold on;
plot(q_final_hat(:,2),'xr');
title('q_1');
R2q2=rsquare(q_final_hat(:,2),q_final(:,2));

subplot 313
plot(q_final(:,3),'o');
hold on;
plot(q_final_hat(:,3),'xr');
title('q_3');
R2q3=rsquare(q_final_hat(:,3),q_final(:,3));


figure;
subplot 311
plot(q_final(:,1),'o');
hold on;
plot(x_final_hat(:,1),'xr');
title('Glucose Concentration');
xlabel('Batch');
ylabel('g/L');
legend('Actual','Predicted');
R2q1=rsquare(x_final_hat(:,1),q_final(:,1));

subplot 312
plot(q_final(:,2),'o');
hold on;
plot(x_final_hat(:,3),'xr');
title('Boimass Concentration');
xlabel('Batch');
ylabel('g/L');
R2q2=rsquare(x_final_hat(:,3),q_final(:,2));

subplot 313
plot(q_final(:,3),'o');
hold on;
plot(x_final_hat(:,4),'xr');
title('Penicillin Concentration');
xlabel('Batch');
ylabel('g/L');
R2q3=rsquare(x_final_hat(:,4),q_final(:,3));

SSM.A = A_hat;
SSM.B = B_hat;
SSM.C = C_hat;
SSM.D = D_hat;
SSM.J_k = J_k;
SSM.PastHorizon = PastHorizon;
SSM.FutureHorizon = FutureHorizon;
SSM.k = k;
SSM.x_typ_idx = x_typ_idx;
SSM.y_typ_idx = y_typ_idx;
SSM.u_typ_idx = u_typ_idx;
SSM.Coeff_Quality = Coeff_Quality;
SSM.mt_final_mean = mt_final_mean;
SSM.mt_final_std = mt_final_std;
SSM.q_final_mean = q_final_mean;
SSM.q_final_std = q_final_std;
SSM.ScaldQualityModel = ScaldQualityModel;


save('SSM_Quality.mat','SSM');