clear
clc;
close all;


load test1;

PastHorizon = 1;
FutureHorizon = 2;

% number of test batches
NTesBatch = 20;

% dimension of the state
k = 7;

% Select the x, y and u to be used in the model
x_typ_idx = [];
y_typ_idx = [2 3 4 5];
u_typ_idx = [1 2 3];

[ P, F, Q, Xu, Yu, Uu ] = UnfoldBatch( bdb.data.x(:,:,1:end-NTesBatch),...
    bdb.data.y(y_typ_idx,:,1:end-NTesBatch), bdb.data.u(u_typ_idx,:,1:end-NTesBatch), PastHorizon, FutureHorizon );

[Pscld, Pmean, Psigma] = zscore(P);
[Fscld, Fmean, Fsigma] = zscore(F);
[Qscld, Qmean, Qsigma] = zscore(Q);
[Yuscld, Ymean, Ysigma] = zscore(Yu);
[Uuscld, Umean, Usigma] = zscore(Uu);



% Pscld1 = ScaleData(P,Pmean, Psigma);
% P1 = ScaleBack(Pscld1, Pmean, Psigma);


[ J, L, D, J_k, SigmaXX, SigmaYY, Omega ] = CVA( Pscld, Fscld, Qscld, k );

[ A_hat, B_hat, C_hat, D_hat ] = SysID( Pscld, Fscld, Qscld, k, Yuscld, Uuscld  );




% select the test batch
i = 50;

% unfold the test data
[ Pt, Ft, Qt, Xt, Yt, Ut ] = UnfoldBatch( bdb.data.x(:,:,i:i),...
    bdb.data.y(y_typ_idx,:,i:i), bdb.data.u(u_typ_idx,:,i:i), PastHorizon, FutureHorizon );

Ptscld = ScaleData(Pt,Pmean,Psigma);

% Y_hatID = Pt*J_k'*C_hat'; 

Y_hatID = Ptscld*J_k'*C_hat'; 

Y_hatID = ScaleBack(Y_hatID,Ymean,Ysigma);


figure;
subplot 221
plot(Yt(1:end,1));
hold on
plot(Y_hatID(1:end,1), 'r--');
title('y_2');
legend('Actual Y','Estimated Y');
R2y1=rsquare(Y_hatID(:,1),Yt(:,1));

subplot 222
plot(Yt(1:end,2));
hold on
plot(Y_hatID(1:end,2), 'r--');
title('y_3');
R2y2=rsquare(Y_hatID(:,2),Yt(:,2));

subplot 223
plot(Yt(1:end,3));
hold on
plot(Y_hatID(1:end,3), 'r--');
title('y_4');
R2y3=rsquare(Y_hatID(:,3),Yt(:,3));

subplot 224
plot(Yt(1:end,4));
hold on
plot(Y_hatID(1:end,4), 'r--');
title('y_5');
R2y4=rsquare(Y_hatID(:,4),Yt(:,4));

% subplot 235
% plot(Yt(1:end,5));
% hold on
% plot(Y_hatID(1:end,5), 'r--');
% title('y_5');
% R2y5=rsquare(Y_hatID(:,5),Yt(:,5));

% subplot 236
% plot(Yt(1:end,6));
% hold on
% plot(Y_hatID(1:end,6), 'r--');
% title('y_6');
% R2y6=rsquare(Y_hatID(:,6),Yt(:,6));

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

save('SSM.mat','SSM');