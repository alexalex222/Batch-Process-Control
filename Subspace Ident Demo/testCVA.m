clear
clc;
close all;

% Set batch system options
opts = BatchSysOpts('nx', 3, 'nu', 3, 'ny', 2, 'inputVar', 'PRBS');

% Generate a system with data
[bdb, LinearSystem] = BatchDataGen(opts, 107, false);

% load test9;

PastHorizon = 1;
FutureHorizon = 2;

% number of test batches
NTesBatch = 3;


[ P, F, Q, Xu, Yu, Uu ] = UnfoldBatch( bdb.data.x(:,:,1:end-NTesBatch),...
    bdb.data.y(:,:,1:end-NTesBatch), bdb.data.u(:,:,1:end-NTesBatch), PastHorizon, FutureHorizon );

% Beta_ARX = ols([P,Q], F);
[~,~,~,~,Beta_ARX] = plsregress([P,Q], F, size([P,Q],2));

Omega = Beta_ARX(end-size(Q,2)+1:end,:);

X=P;
Y=F;
k=3;

% [~,~,XS,~,Beta_PLS]=plsregress(P,F,15);


[r,p]=corrcoef(zscore(P));
[i,j] = find(p<0.05);

[sample,dx] = size(X);
[~,dy] = size(Y);

% SigmaXX = cov(X);
SigmaXX = (X-repmat(mean(X),sample,1))'*(X-repmat(mean(X),sample,1))/(sample-1);
SigmaYY = cov(Y);
SigmaXY = (X-repmat(mean(X),sample,1))'*(Y-repmat(mean(Y),sample,1))/(sample-1);

SigmaXX_o=SigmaXX;

SigmaXX = nearestSPD(SigmaXX);



[U, D, V] = svd(SigmaXX^(-0.5)*SigmaXY*SigmaYY^(-0.5));

J = U'*SigmaXX^(-0.5);
J_k = U(:,1:k)'*SigmaXX^(-0.5);
L = V'*SigmaYY^(-0.5);





B = J*J'*(X-repmat(mean(X),sample,1))'*(Y-repmat(mean(Y),sample,1));
F_hat = (X-repmat(mean(X),sample,1))*B+repmat(mean(Y),sample,1);


