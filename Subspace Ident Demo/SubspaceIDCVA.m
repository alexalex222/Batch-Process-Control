clear
clc;
close all;

% Set batch system options
opts = BatchSysOpts('nx', 3, 'nu', 3, 'ny', 2, 'inputVar', 'PRBS',...
    'inputTrend','randomwalk','meas_STD',0.3);

% Generate a system with data
[bdb, LinearSystem] = BatchDataGen(opts, 106, true);

% select the variables to be used
x_typ_idx = [];
y_typ_idx = [];
u_typ_idx = [];

% Determine the past and future horizon
PastHorizon = 2;
FutureHorizon = 2;

% number of test batches
NTesBatch = 1;
NTrainBatch = 30;

[ P, F, Q, Xu, Yu, Uu, nb ] = UnfoldBatch( bdb.data.x(:,:,1:NTrainBatch),...
    bdb.data.y(:,:,1:NTrainBatch), bdb.data.u(:,:,1:NTrainBatch), PastHorizon, FutureHorizon );

% the dimension of the estimated state
k = 3;

[ J, L, D, J_k, SigmaXX, SigmaYY ] = CVA( P, F, Q, k );

[ A_hat, B_hat, C_hat, D_hat ] = SysID( P, F, Q, k, Yu, Uu  );

[sample,dp] = size(P);


mt = P*J_k';


figure;

subplot 311
plot(Xu(1:500,1));
hold on
plot(mt(1:500,1), 'r--');
title('x_1');
I = legend('$X$','$\hat{X}$');
set(I,'Interpreter','Latex');

subplot 312
plot(Xu(1:500,2));
hold on
plot(mt(1:500,2), 'r--');
title('x_2');
% legend('Estimated State','Actual State');



subplot 313
plot(Xu(1:500,3));
hold on
plot(mt(1:500,3), 'r--');
title('x_3');
% legend('Estimated State','Actual State');






% unfold the test data
[ Pt, Ft, Qt, Xt, Yt, Ut ] = UnfoldBatch(  bdb.data.x(:,:,end-NTesBatch+1:end),...
    bdb.data.y(:,:,end-NTesBatch+1:end), bdb.data.u(:,:,end-NTesBatch+1:end), PastHorizon, FutureHorizon );

mt_t = Pt*J_k';

% Y_hatID = (Pt-repmat(mean(P),size(Pt,1),1))*J_k'*C_hat'+Ut*D_hat';
Y_hatID = Pt*J_k'*C_hat'+Ut*D_hat';

figure;
subplot 211
plot(Yt(1:end,1));
hold on
plot(Y_hatID(1:end,1), 'r--');
title('y_1');
legend('Actual Y','Estimated Y');
R2y1=rsquare(Y_hatID(:,1),Yt(:,1));

subplot 212
plot(Yt(1:end,2));
hold on
plot(Y_hatID(1:end,2), 'r--');
title('y_2');
legend('Actual Y','Estimated Y');
R2y2=rsquare(Y_hatID(:,2),Yt(:,2));


% find the relationship between the actual states and identified states
% [~,~,~,~,Beta_PLS]=plsregress(Xu,mt,3);
% X_hatpls=[ones(size(Xu,1),1),Xu]*Beta_PLS;
% figure;
% subplot 311
% plot(X_hatpls(1:300,1), 'r--');
% hold on
% plot(mt(1:300,1));
% title('x_1');
% legend('Estimated Y','Actual Y');
% R2m1=rsquare(X_hatpls(:,1),mt(:,1));
% 
% subplot 312
% plot(X_hatpls(1:300,2), 'r--');
% hold on
% plot(mt(1:300,2));
% title('x_2');
% legend('Estimated Y','Actual Y');
% R2m2 = rsquare(X_hatpls(:,2),mt(:,2));
% 
% subplot 313
% plot(X_hatpls(1:300,3), 'r--');
% hold on
% plot(mt(1:300,3));
% title('x_3');
% legend('Estimated Y','Actual Y');
% R2m3 = rsquare(X_hatpls(:,3),mt(:,3));

% [~,~,~,~,Beta_PLS]=plsregress(mt,Xu,3);
% X_hatpls=[ones(size(mt,1),1),mt]*Beta_PLS;

P1 = ols(Xu, mt);
% P2 = ols(mt, Xu);
P2 = ols(mt(96:96:end,:), Xu(96:96:end,:));


X_hatpls = mt_t*P2;


figure;
subplot 311
plot(Xt(1:end,1));
hold on
plot(X_hatpls(1:end,1), 'r--');

title('x_1');
I = legend('$X$','$\hat{X}R$');
set(I,'Interpreter','Latex');

R2m1=rsquare(X_hatpls(:,1),Xt(:,1));

subplot 312
plot(Xt(1:end,2));
hold on
plot(X_hatpls(1:end,2), 'r--');
title('x_2');
% legend('Estimated X','Actual X');
R2m2 = rsquare(X_hatpls(:,2),Xt(:,2));

subplot 313
plot(Xt(1:end,3));
hold on
plot(X_hatpls(1:end,3), 'r--');
title('x_3');
R2m3 = rsquare(X_hatpls(:,3),Xt(:,3));



% figure;
% subplot 311
% plot(Xu(1:300,1));
% hold on
% plot(X_hatpls(1:300,1), 'r--');
% 
% title('x_1');
% I = legend('$X$','$\hat{X}R$');
% set(I,'Interpreter','Latex');
% 
% R2m1=rsquare(X_hatpls(:,1),Xu(:,1));
% 
% subplot 312
% plot(Xu(1:300,2));
% hold on
% plot(X_hatpls(1:300,2), 'r--');
% title('x_2');
% % legend('Estimated X','Actual X');
% R2m2 = rsquare(X_hatpls(:,2),Xu(:,2));
% 
% subplot 313
% plot(Xu(1:300,3));
% hold on
% plot(X_hatpls(1:300,3), 'r--');
% 
% title('x_3');
% legend('Estimated X','Actual X');
% R2m3 = rsquare(X_hatpls(:,3),Xu(:,3));


% figure;
% subplot 311
% plot(Uu(1:500,1));
% title('u_1');
% % legend('Estimated Y','Actual Y');
% 
% subplot 312
% plot(Uu(1:500,2));
% title('u_2');
% % legend('Estimated Y','Actual Y');
% 
% subplot 313
% plot(Uu(1:500,3));
% title('u_3');

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
save('bdb.mat','bdb');

%% PLS regression model
% [~,~,~,~,Beta_PLS]=plsregress(P,F(:,1:2),3);
% Y_hatpls=[ones(size(Pt,1),1),Pt]*Beta_PLS;
% 
% figure;
% subplot 211
% plot(Y_hatpls(1:end,1), 'r--');
% hold on
% plot(Yt(1:end,1));
% title('y_1');
% legend('Estimated Y','Actual Y');
% R2y11=rsquare(Y_hatpls(:,1),Yt(:,1));
% 
% subplot 212
% plot(Y_hatpls(1:end,2), 'r--');
% hold on
% plot(Yt(1:end,2));
% title('y_2');
% legend('Estimated Y','Actual Y');
% R2y22=rsquare(Y_hatpls(:,2),Yt(:,2));
