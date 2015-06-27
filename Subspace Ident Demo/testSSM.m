clear;
clc;
close all;

load bdb;
load SSM;

A = SSM.A;
B = SSM.B;
C = SSM.C;
D = SSM.D;
J_k = SSM.J_k;

PastHorizon = 2;
FutureHorizon = 2;
NTesBatch = 1;

% select a test batch
i = 60;

% X = bdb.data.x(2:5,:,i)';
% Y = bdb.data.y(:,:,i)';
% U = bdb.data.u(1:3,:,i)';

[ P, F, Q, Xu, Yu, Uu ] = UnfoldBatch( bdb.data.x(:,:,1:end-NTesBatch),...
    bdb.data.y(:,:,1:end-NTesBatch), bdb.data.u(:,:,1:end-NTesBatch), PastHorizon, FutureHorizon );

[ Pt, Ft, Qt, Xt, Yt, Ut ] = UnfoldBatch( bdb.data.x(:,:,i:i),...
    bdb.data.y(:,:,i:i), bdb.data.u(:,:,i:i), PastHorizon, FutureHorizon );

Y = bdb.data.y(:,:,i);
U = bdb.data.u(:,:,i);
U(:,end+1) = U(:,end);
Y_hat =zeros(size(Y));

% set the current sampling instance
k = 3;
Y_hat(:,1:k-1) = Y(:,1:k-1);

for j = k:101
    pt = [Y_hat(:,j-2)' Y_hat(:,j-1)' U(:,j-2)' U(:,j-1)'];
%     pt = [ Y_hat(:,j-1)'  U(:,j-1)'];
    Y_hat(:,j) = SSM.C*SSM.J_k*pt';
end

Coeff = SSM.C*SSM.J_k;

Y = Y';
Y_hat = Y_hat';

figure;
subplot 211
plot(Y(:,1));
hold on
plot(Y_hat(:,1), 'r--');
title('y_1');
legend('Actual Y','Estimated Y');
R2y1=rsquare(Y_hat(:,1),Y(:,1));

subplot 212
plot(Y(:,2));
hold on
plot(Y_hat(:,2), 'r--');
title('y_2');
R2y2=rsquare(Y_hat(:,2),Y(:,2));

% subplot 223
% plot(Y(:,3));
% hold on
% plot(Y_hat(:,3), 'r--');
% title('y_3');
% R2y3=rsquare(Y_hat(:,3),Y(:,3));
% 
% subplot 224
% plot(Y(:,4));
% hold on
% plot(Y_hat(:,4), 'r--');
% title('y_4');
% R2y4=rsquare(Y_hat(:,4),Y(:,4));

