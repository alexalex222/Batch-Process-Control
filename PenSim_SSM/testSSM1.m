clear;
clc;
close all;

load test9;
load SSM;



A = SSM.A;
B = SSM.B;
C = SSM.C;
D = SSM.D;
J_k = SSM.J_k;

PastHorizon = 1;
FutureHorizon = 2;
NTesBatch = 1;


% Select the x, y and u to be used in the model
x_typ_idx = SSM.x_typ_idx;
y_typ_idx = SSM.y_typ_idx;
u_typ_idx = SSM.u_typ_idx;

[ P, F, Q, Xu, Yu, Uu ] = UnfoldBatch( bdb.data.x(:,:,1:end-NTesBatch),...
    bdb.data.y(y_typ_idx,:,1:end-NTesBatch), bdb.data.u(u_typ_idx,:,1:end-NTesBatch), PastHorizon, FutureHorizon );

[Pscld, Pmean, Psigma] = zscore(P);
[Fscld, Fmean, Fsigma] = zscore(F);
[Qscld, Qmean, Qsigma] = zscore(Q);
[Yuscld, Ymean, Ysigma] = zscore(Yu);
[Uuscld, Umean, Usigma] = zscore(Uu);

Coeff_ARX = ols(Pscld,Fscld);


% select a test batch
i = 60;

Y = bdb.data.y(2:5,:,i)';
U = bdb.data.u(1:3,:,i)';
U(end+1,:) = U(end,:)';
Y_hat =zeros(size(Y));
Y_hat_scld =zeros(size(Y));

Y_scld = ScaleData(Y, Ymean, Ysigma);
U_scld = ScaleData(U, Umean, Usigma);

% set the current sampling instance
k = 3;
Y_hat(1:k-1,:) = Y(1:k-1,:);
Y_hat_scld(1:k-1,:) = Y_scld(1:k-1,:);

for j = k:601
%     pt = [Y_hat_scld(j-2,:) Y_hat_scld(j-1,:) U_scld(j-2,:) U_scld(j-1,:)];
    pt = [ Y_hat_scld(j-1,:)  U_scld(j-1,:)];
    temp = (SSM.C*SSM.J_k*pt')';
    
    Y_hat_scld(j,:) = temp;
    Y_hat(j,:) = ScaleBack(temp,Ymean, Ysigma);
    
end

Coeff = SSM.C*SSM.J_k;

% Y = Y';
% Y_hat = Y_hat';

figure;
subplot 221
plot(Y(:,1));
hold on
plot(Y_hat(:,1), 'r--');
title('y_1');
legend('Actual Y','Estimated Y');
R2y1=rsquare(Y_hat(:,1),Y(:,1));

subplot 222
plot(Y(:,2));
hold on
plot(Y_hat(:,2), 'r--');
title('y_2');
R2y2=rsquare(Y_hat(:,2),Y(:,2));

subplot 223
plot(Y(:,3));
hold on
plot(Y_hat(:,3), 'r--');
title('y_3');
R2y3=rsquare(Y_hat(:,3),Y(:,3));

subplot 224
plot(Y(:,4));
hold on
plot(Y_hat(:,4), 'r--');
title('y_4');
R2y4=rsquare(Y_hat(:,4),Y(:,4));

