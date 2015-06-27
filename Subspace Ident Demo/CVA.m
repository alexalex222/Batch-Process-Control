function [ J, L, D, J_k, SigmaXX, SigmaYY, Omega ] = CVA( X, Y, Q, k )
%Canonical Variate Analysis
%   Input
%   X: past input/outputs
%   Y: future outputs
%   Q; future inputs
%   k: the dimension of subspace
%   Output
%   J: loading matrix for X
%   L: loading matrix for Y
%   D: diagonal matrix containing canonical correaltion

[sample,dx] = size(X);
% [~,dy] = size(Y);


% Future data correction
% Beta_ARX = ols([X,Q], Y);

[~,~,~,~,Beta_ARX] = plsregress([X,Q], Y, size([X,Q],2));
Omega = Beta_ARX(end-size(Q,2)+1:end,:);
Y = Y - Q*Omega;

% SigmaXX = cov(X);
SigmaXX = (X-repmat(mean(X),sample,1))'*(X-repmat(mean(X),sample,1))/(sample-1);
SigmaYY = cov(Y);
SigmaXY = (X-repmat(mean(X),sample,1))'*(Y-repmat(mean(Y),sample,1))/(sample-1);

SigmaXX = nearestSPD(SigmaXX);
SigmaYY = nearestSPD(SigmaYY);

[U, D, V] = svd(SigmaXX^(-0.5)*SigmaXY*SigmaYY^(-0.5));

J = U'*SigmaXX^(-0.5);
J_k = U(:,1:k)'*SigmaXX^(-0.5);
L = V'*SigmaYY^(-0.5);


end

